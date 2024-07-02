'''
Created on 28.11.2016

@author: witek
'''

import itertools
import copy

class Graph:
    def __init__(self,edges):
        neighbourhoods=[]
        neighbourhoods.extend(edges)
        neighbourhoods.extend((y,x) for (x,y) in edges)
        
        neighbourhoods.sort()
       
        self.lists={}
        for k,g in itertools.groupby(neighbourhoods,key=lambda x:x[0]):
            self.lists[k]=list(y for (x,y) in g)
        
    def degeneracy_order(self):
        tmp=copy.deepcopy(self.lists)
              
        result=[]
        
        while len(tmp)>0:
            
            (node,_)=min(((k,len(tmp[k])) for k in tmp.keys()),key=lambda x:x[1])            
            
            tmp.pop(node)
            
            for x in tmp.keys():
                try:
                    tmp[x].remove(node)
                except ValueError:
                    pass

            result.append(node)
            
        return result
      
    def __bron_kerbosch_rek(self,R,P,X,result):
        if len(P)==0 and len(X)==0:
            result.append(frozenset(R))
            return
        
        pivot_candidates=P.union(X)
        (pivot_node,_)=max(((node,len(self.lists[node])) for node in pivot_candidates),key=lambda x:x[1])
        
        for node in P.difference(self.lists[pivot_node]):
            self.__bron_kerbosch_rek(R.union({node}),P.intersection(self.lists[node]),X.intersection(self.lists[node]),result)
            P.difference_update({node})
            X.update({node})
    
    def bron_kerbosch(self):
        
        order=self.degeneracy_order()
        P=set(self.lists.keys())
        R=set()
        X=set() 
        result=[]
                
        for node in order:
            self.__bron_kerbosch_rek(R.union({node}),P.intersection(self.lists[node]),X.intersection(self.lists[node]),result)
            P.difference_update({node})
            X.update({node})
        
        return result