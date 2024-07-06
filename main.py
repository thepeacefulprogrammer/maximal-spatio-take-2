'''
Created on 28.11.2016

@author: witek
'''

import structures.input
from algorithms.mainalg import maxspatiotempcolloc
from algorithms.bb_randy_impl import BBmaxspatiotempcolloc
from algorithms.mdcop import mdcop
from algorithms.mdcop import find_maximal
from algorithms.utils import manhattan_distance
from time import perf_counter
from experiments.batchgenerator import I,C,P,executeAll
import csv
import os
from multiprocessing import Process, Queue
import psutil
import logging
import platform
from structures.graph import Graph


# Override the build_local_candidates method
def build_local_candidates_override(self):
    g = Graph(self.prevalent_pairs)
    self.local_candidates = g.bron_kerbosch()

    if self.local_candidates:
        lenc = max(len(c) for c in self.local_candidates)
        self.local_candidates_divided = [None for _ in range(lenc + 1)]
        while lenc >= 2:
            self.local_candidates_divided[lenc] = set(c for c in self.local_candidates if len(c) == lenc)
            lenc = lenc - 1
    else:
        self.local_candidates_divided = []
    
    logging.info(f"Built local candidates. Count: {len(self.local_candidates)}")

# Apply the override
structures.input.WorldState.build_local_candidates = build_local_candidates_override


if platform.system() != 'Windows':
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        print(f"File descriptor limits - Soft: {soft}, Hard: {hard}")
    except ImportError:
        print("The 'resource' module is not available. Cannot get file descriptor limits.")
else:
    print("Running on Windows. File descriptor limits are not applicable.")



def basic_test():
    data=structures.input.InputDataset()
    data.loadFromFile("large_testMDCOP.txt")
    
    new_alg_result=maxspatiotempcolloc(data,1.5,0.5,2/3,False)    
    print("New algorithm result (original): ",new_alg_result)
        
    old_alg_result=mdcop(data,1.5,0.5,2/3,False)
    print("Old algorithm result (non maximal): ",old_alg_result)
    
    converted_result=sorted(convert_set_to_tuple(new_alg_result))
    print("New algorithm result (common form): ",converted_result)
    
    max_old_alg_result=sorted(find_maximal(old_alg_result))
    print("Old algorithm result (common form): ",max_old_alg_result)
    

class CSVResultSaver:
    def __init__(self,file):
        self.file=file
        self.writer = csv.writer(file, delimiter=';', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        
    def feed(self,data):
        row=list(data[0:2])
        row.extend(list(data[2]))
        self.writer.writerow(row)
        self.file.flush()
        os.fsync(self.file.fileno())        


def experiment_all(filename,minprev,minfreq,maxdist,verbose):
    print("Experimenting...")
    
    data=structures.input.InputDataset()
    data.loadFromFile(filename)
    t1=perf_counter()    
    res1=maxspatiotempcolloc(data,maxdist,minprev,minfreq,False,False,verbose)
    t2=perf_counter()
    print("New algorithm, no predictions, no tree saving: ",t2-t1)
    
    data=structures.input.InputDataset()
    data.loadFromFile(filename)
    t3=perf_counter()
    res2=maxspatiotempcolloc(data,maxdist,minprev,minfreq,True,False,verbose)    
    t4=perf_counter()
    print("New algorithm predictions, no tree saving: ",t4-t3)
    
    data=structures.input.InputDataset()
    data.loadFromFile(filename)
    t5=perf_counter()
    res3=maxspatiotempcolloc(data,maxdist,minprev,minfreq,False,True,verbose)
    t6=perf_counter()
    print("New algorithm, no predictions, tree saving: ",t6-t5)
    
    data=structures.input.InputDataset()
    data.loadFromFile(filename)
    t7=perf_counter()
    res4=maxspatiotempcolloc(data,maxdist,minprev,minfreq,True,True,verbose)    
    t8=perf_counter()
    print("New algorithm predictions, tree saving: ",t8-t7)
    
    data=structures.input.InputDataset()
    data.loadFromFile(filename)
    t9=perf_counter()
    res5=maxspatiotempcolloc(data,maxdist,minprev,minfreq,False,True,verbose)
    t10=perf_counter()
    print("New algorithm, no predictions, tree saving, tree cleaning: ",t10-t9)
    
    data=structures.input.InputDataset()
    data.loadFromFile(filename)
    t11=perf_counter()
    res6=maxspatiotempcolloc(data,maxdist,minprev,minfreq,True,True,verbose)    
    t12=perf_counter()
    print("New algorithm predictions, tree saving, tree cleaning: ",t12-t11)
    
    data=structures.input.InputDataset()
    data.loadFromFile(filename)    
    t13=perf_counter()
    res7=mdcop(data,maxdist,minprev,minfreq,verbose)
    t14=perf_counter()
    print("Old algorithm: ",t14-t13)  

    data = structures.input.InputDataset()
    data.loadFromFile(filename)
    t15 = perf_counter()
    res8 = BBmaxspatiotempcolloc(data, maxdist, minprev, minfreq, manhattan_distance, False, verbose, False)
    t16 = perf_counter()
    print("New bounding box algorithm: ", t16 - t15)

    print("Verifying...")
    res1a=sorted(convert_set_to_tuple(res1))
    res2a=sorted(convert_set_to_tuple(res2))
    res3a=sorted(convert_set_to_tuple(res3))
    res4a=sorted(convert_set_to_tuple(res4))
    res5a=sorted(convert_set_to_tuple(res5))
    res6a=sorted(convert_set_to_tuple(res6))
    res7a=sorted(find_maximal(res7))
    comOldvsNew1 = (res7a==res1a)
    comOldvsNew2 = (res7a==res2a)
    comOldvsNew3 = (res7a==res3a)
    comOldvsNew4 = (res7a==res4a)
    comOldvsNew5 = (res7a==res5a)
    comOldvsNew6 = (res7a==res6a)
    print("Results equal:")
    print(" Old vs no predictions, no tree saving: ", comOldvsNew1)
    print(" Old vs predictions, no tree saving: ", comOldvsNew2)
    print(" Old vs no predictions, tree saving: ", comOldvsNew3)
    print(" Old vs predictions, tree saving: ", comOldvsNew4)
    print(" Old vs no predictions, tree saving, tree cleaning: ", comOldvsNew5)
    print(" Old vs predictions, tree saving, tree cleaning: ", comOldvsNew6)
    
    return filename,maxdist,minprev,minfreq,t2-t1,t4-t3,t6-t5,t8-t7,t10-t9,t12-t11,t14-t13,comOldvsNew1, comOldvsNew2,comOldvsNew3, comOldvsNew4,comOldvsNew5, comOldvsNew6
def get_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    # Return memory usage in MB
    return mem_info.rss / (1024 * 1024)


def experiment_one_subprocess(q, filename, minprev, minfreq, maxdist, verbose, alg_name):
    mem = psutil.virtual_memory()

    data = structures.input.InputDataset()
    data.loadFromFile(filename)

    print("Experimenting on ", alg_name)
    t1 = perf_counter()
    res = []  # Initialize res with a default value
    if alg_name == "MAXMDCOP_PRED=0_SAVE=0":
        res = maxspatiotempcolloc(data, maxdist, minprev, minfreq, False, False, verbose, False)
    elif alg_name == "MAXMDCOP_PRED=1_SAVE=0":
        res = maxspatiotempcolloc(data, maxdist, minprev, minfreq, True, False, verbose, False)
    elif alg_name == "MAXMDCOP_PRED=0_SAVE=1":
        res = maxspatiotempcolloc(data, maxdist, minprev, minfreq, False, True, verbose, False)
    elif alg_name == "MAXMDCOP_PRED=1_SAVE=1":
        res = maxspatiotempcolloc(data, maxdist, minprev, minfreq, True, True, verbose, False)
    elif alg_name == "MAXMDCOP_PRED=0_SAVE=2":
        res = maxspatiotempcolloc(data, maxdist, minprev, minfreq, False, True, verbose, True)
    elif alg_name == "MAXMDCOP_PRED=1_SAVE=2":
        res = maxspatiotempcolloc(data, maxdist, minprev, minfreq, True, True, verbose, True)
    elif alg_name == "FASTMDCOPMINER":
        res = mdcop(data, maxdist, minprev, minfreq, verbose)
    elif alg_name == "BB_RANDY":
        res = BBmaxspatiotempcolloc(data, maxdist, minprev, minfreq, manhattan_distance, False, verbose, False)
    else:
        print("Unknown algorithm")
    t2 = perf_counter()

    peak_memory_usage = get_memory_usage()

    print("Measured time: ", t2 - t1)
    print("Available memory at start: ", mem.available / (1024 * 1024), "MB")
    print("Peak memory usage: ", peak_memory_usage, "MB")

    q.put((filename, alg_name, maxdist, minprev, minfreq, t2 - t1, mem.available / (1024 * 1024), peak_memory_usage, res))

def experiment_one(filename,minprev,minfreq,maxdist,verbose,alg_name):
    queue = Queue()
    p = Process(target=experiment_one_subprocess, args=(queue,filename,minprev,minfreq,maxdist,verbose,alg_name))
    p.start()
    p.join() # this blocks until the process terminates
    return queue.get()



experiment_all_tasks= [
    {
        'filename':C('large_testMDCOP.txt'),
        'maxdist':I([5]),
        'minfreq':I([0.9]),
        'minprev':I([0.5, 0.7]),
        'verbose':C(0),
    }
    ]


def perform_experiment_all(output_name,tasks):
    with open(output_name, 'w') as csvfile:
        result=executeAll(tasks,experiment_all, False, CSVResultSaver(csvfile))
        
experiment_one_tasks=[
    {
        'filename':C('pigeons_zonnani_rd2_2880T.csv'),
        #'filename':C('large_testMDCOP.txt'),
        #'filename':C('testMDCOP21.txt'),
        # 'filename':C('synthesized_data.csv'),
        'maxdist':I([5]),
        'minfreq':I([0.9]),
        'minprev':I([0.5, 0.7]),
        'verbose':C(1),
        'alg_name':I(["BB_RANDY","MAXMDCOP_PRED=0_SAVE=0","MAXMDCOP_PRED=1_SAVE=0","MAXMDCOP_PRED=0_SAVE=1","MAXMDCOP_PRED=1_SAVE=1","MAXMDCOP_PRED=0_SAVE=2","MAXMDCOP_PRED=1_SAVE=2","FASTMDCOPMINER"])

    }
    ]

import random

def generate_patterned_dataset(filename, num_records, max_timestamps=10, max_object_ids=6):
    with open(filename, 'w') as file:
        for i in range(num_records):
            timestamp = (i // max_object_ids) + 1  # Sequential timestamps
            obj_id = (i % max_object_ids) + 1  # Sequential object IDs, limited to max_object_ids
            flag = 1  # Always 1

            # Following the identified pattern
            val1 = (i // max_object_ids) % 7  # Value 1 based on timestamp
            val2 = (i % max_object_ids) % 7  # Value 2 based on object ID

            file.write(f"{timestamp},{obj_id},{flag},{val1},{val2}\n")



def perform_experiment_one(output_name,tasks):
    with open(output_name, 'w') as csvfile:
        result=executeAll(tasks,experiment_one, False, CSVResultSaver(csvfile))       

if __name__ == '__main__':       
        #perform_experiment_all("result5.csv",experiment_all_tasks)
        #generate_patterned_dataset('large_testMDCOP.txt', 500000)
        perform_experiment_one("result4.csv",experiment_one_tasks)
