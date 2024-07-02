'''
Created on 03.01.2017

@author: witek
'''
import inspect
import copy

def I(data):
    return (1,data)

def P(data):
    return (2,data)

def C(data):
    return (3,data)

def expandIterables(iterables,pos,data,out):
    if pos<len(iterables):
        name=iterables[pos][0]
        iterable=iterables[pos][1]
        for x in iterable:
            data[name]=x
            expandIterables(iterables, pos+1, data, out)
    else:
        out.append(copy.deepcopy(data))

def expand(tasks):
    result=[]    
    for block in tasks:
        iterables=[]
        procedures=[]
        data={}
            
        for key, value in block.items():                                    
            if (value[0]==1):
                iterables.append((key,value[1]))
            if (value[0]==2):
                procedures.append((key,value[1]))
            if (value[0]==3):
                data[key]=value[1]
                
        block_result=[]        
        expandIterables(iterables, 0, data, block_result)
        
        for x in procedures:
            args=inspect.getargspec(x[1]).args
            call={}
            
            for data in block_result:
                for arg in args:
                    call[arg]=data[arg]
                
                data[x[0]]=x[1](**call)

        
        result.append(block_result)
        
    return result


def executeAll(tasks,fun,mapped_args=True,result_saver=None):
    results=[]
    expanded_tasks=expand(tasks)
    for blocknum in range(0,len(expanded_tasks)):
        print("Block: ",blocknum+1,"/",len(expanded_tasks))
        subblock=expanded_tasks[blocknum]
        for subblocknum in range(0,len(subblock)):
            print("\tSubblock: ",subblocknum+1,"/",len(subblock))
            data=subblock[subblocknum]
            for x in sorted(data.keys()):
                print("\t\t",x,"=",data[x])
            
            if fun is not None:
                if mapped_args:
                    funres=fun(data)
                else:
                    funres=fun(**data)
                    
                res=(blocknum+1,subblocknum+1,funres)
                results.append(res)
                                    
                if result_saver is not None:
                        result_saver.feed(res)
    
    return results


'''
tasks=[
        {        
        'a':I([1,2,3]),
        'b':P((lambda a:'a'+str(a)+'x')),
        'c':C('test')
        },       
        {        
        'a':I([1,2,3]),
        'e':P((lambda a:'y'+str(a)+'y')),
        'f':C(666)
        }
        
   ]


#print(expand(tasks))
#executeAll(tasks,None)


tasks2=[
    {
        'x':I([1,2,3]),
        'y':I([4,5,6]),
        'z':P(lambda x,y:x+y)
    },
    {
        'x':I([1,2,3]),
        'y':I([4,5,6]),
        'z':P(lambda x,y:x-y)
    }
    ]


def fun1(x,y,z):
    print(x,y,z)
    return 1

print(executeAll(tasks2,fun1,False))
'''
