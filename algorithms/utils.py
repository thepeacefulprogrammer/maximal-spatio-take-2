'''
Created on 29.11.2016

@author: witek
'''

import itertools
import functools

'''
Simple map reduce function. Executes map_func and key_func on input data from iterable. Performs reduction w.r.t the computed key using reduce_func.
'''
def map_reduce(iterable,key_func,map_func,reduce_func):
    data=[(key_func(x),map_func(x)) for x in iterable]
    data.sort(key=lambda x:x[0], reverse=False)
    
    print(data)
    
    result={}
    
    for k,g in itertools.groupby(data,key=lambda x:x[0]):                
        result[k]=functools.reduce(lambda x,y:reduce_func(x,y[1]), g,next(g)[1])
        
    return result
        
'''
Empty iterator - returns nothing
'''
def empty_iterator():
    return
    yield


class PartList:
    '''
    Iterator over part of a list - from index start iterate over count elements. If list ends sooner then iterator finishes as well. 
    MAYBE its faster than [x:y]. Most probably its the same or slower than generator expression. Candidate for removal.
    '''
    
    def __init__(self,list,start,count):
        self.list=list
        self.start=start
        self.end=start+count
        
    def __iter__(self):
        return self
    
    def __next__(self):
        if self.start>=min(len(self.list),self.end):
            raise StopIteration
        
        self.start+=1
        
        return self.list[self.start-1]
    


    
