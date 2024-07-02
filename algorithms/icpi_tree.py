'''
Created on 28.11.2016

@author: witek
'''

import itertools
from algorithms.utils import empty_iterator
from algorithms.utils import PartList
from itertools import groupby

class ICPITree:

    def __init__(self):
        self.neighbors=None
        self.dictionary=None
        self.neighborlist=None
        
    def build(self,neighbors):
        sortedneighbors=sorted(neighbors)#,key=lambda x:(x[0].feature,x[0].instance,x[1].feature,x[1].instance))
        self.neighbors=neighbors
        
        self.dictionary={}
        
        pos=0
        for  k,g in itertools.groupby(sortedneighbors,key=lambda x:(x[0],x[1].feature)):
            count=0
            for _ in g:
                count+=1                            
            self.dictionary[k]=(pos,count)
            pos+=count
        
        self.neighborlist=[x for (y,x) in sortedneighbors]
        
        #print(self.dictionary)
            
    def find_neighbors(self,fid,feature):
        try:
            (pos,count)=self.dictionary[(fid,feature)]
        except:
            return empty_iterator()
        #return PartList(self.neighborlist,pos,count)
        #return self.neighborlist[pos:pos+count]
        return (self.neighborlist[p] for p in range(pos,min(len(self.neighborlist),pos+count)))

    def neighbor_count(self,fid,feature):
        (_,count)=self.dictionary[(fid,feature)]
        return count
    
    
    def prevalent_pairs(self,minprev,instance_counts):
        
        sortedneighbors=sorted(self.neighbors,key=lambda x:(x[0].feature,x[1].feature))                
        
        result=set()
        
        for k,g in itertools.groupby(sortedneighbors,key=lambda x:(x[0].feature,x[1].feature)):
            pair_instances=list(g)
            count_p1=len(set([x for (x,y) in pair_instances]))
            count_p2=len(set([y for (x,y) in pair_instances]))
            prevalence=min(count_p1/instance_counts[k[0]],count_p2/instance_counts[k[1]]);
            if prevalence>=minprev: 
                result.add(k)
        
        return result
    
    def prevalent_pairs_with_instances(self,minprev,instance_counts,forbidden_pairs):
        #tutaj robie podobne rzeczy jak procedura powyżej, ale korzystam z drzewa icpi. Oszczędzam trochę na czasie sortowania
        result=[]
        instances=[]
        
        keys=sorted( self.dictionary.keys(),key=lambda x: (x[0].feature,x[1],x[0].instance))
        for pair,g in groupby(keys,key=lambda x: (x[0].feature,x[1])):
            if pair in forbidden_pairs:
                continue
            
            pair_instances=[]            
            for k in g:
                n1=k[0]
                (pos,count)=self.dictionary[k]
                for n2 in self.neighborlist[pos:pos+count]:
                    pair_instances.append((n1,n2))
                            
            count_p1=len(set([x for (x,y) in pair_instances]))
            count_p2=len(set([y for (x,y) in pair_instances]))
            prevalence=min(count_p1/instance_counts[pair[0]],count_p2/instance_counts[pair[1]]);
            if prevalence>=minprev: 
                result.append(pair)
                instances.append(pair_instances)                
        
        return (result,instances)