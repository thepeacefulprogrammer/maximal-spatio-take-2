'''
Created on 28.11.2016

@author: witek
'''

import math
import csv
from algorithms.icpi_tree import ICPITree
from algorithms.planesweep import planesweep
from structures.graph import Graph
from itertools import filterfalse

#from sortedcontainers import SortedDict
#from bisect import bisect_left

class Position:
    '''
    Represents position in 2D cartesian coordinates
    '''    
    
    def __init__(self,x=0,y=0):
        '''
        Allows to specify position. Default position is 0,0.
        '''
        self.x=x
        self.y=y
        
    def distance(self,other):
        '''
        Computes euclidean distance of the position to another position
        '''
        return math.sqrt(math.pow(self.x-other.x,2)+math.pow(self.y-other.y,2))
    
    

class FeatureInstanceIdentifier:
    '''
    Represents feature instance identifier which is composed of feature number and instance number.
    '''    
    def __init__(self, feature, instance):
        self.feature=feature
        self.instance=instance
    
    def __hash__(self):
        return hash((self.feature,self.instance))
    
    def __eq__(self,other):
        return (self.feature,self.instance)==(other.feature,other.instance)
    
    def __str__(self):
        return "<"+self.feature.__str__()+","+self.instance.__str__()+">";
    
    def __repr__(self):
        return self.__str__()
    
    def __lt__(self,other):
        return ((self.feature,self.instance)<(other.feature,other.instance))
        
    
class FeatureInstance():
    '''
    Represents a feature instance. It has a position as well as feature instance identifier.
    '''
    def __init__(self, position, instanceid):
        self.pos=position
        self.id=instanceid
        
                  
class TrieNode():
    def __init__(self):
        self.children={}
        self.value=None
        
    def __insert__(self,x,p,v):                
        if p<len(x):
            child_node=self.children.get(x[p])
            if child_node is None:
                child_node=TrieNode()
                self.children[x[p]]=child_node
                
            child_node.__insert__(x,p+1,v)
        else:
            self.value=v
                
    def insert(self,x,v):
        self.__insert__(x,0,v)
        
    def __get__(self,x,p):
        if p<len(x):
            child_node=self.children.get(x[p])
            if child_node is None:
                return None
                            
            return child_node.__get__(x,p+1)
        else:
            child_node=self
            while child_node.value is None:                
                child_node=next(iter(child_node.children.values()))
            #zwrocenie wartosci a nie kopii w tym momencie jest teoretycznie błędne, bo 
            #gdyby cała wyszukana sekwencja miała być prefiksem nowej to uszkodzilibyśmy tablicę poziomów tej sekwencji. 
            #Coś takiego jednak nie zajdzie bo nie wydłużamy istniejących całych sekwencji - idziemy "od góry"
            #Na wszelki wypadek jednak zwracam kopię
            return child_node.value[0:len(child_node.value)]
        
    def get(self,x):
        return self.__get__(x,0)
    
    def __get_prefix_length__(self,x,p):
        if p<len(x):
            child_node=self.children.get(x[p])
            if child_node is None:
                return p
                            
            return child_node.__get_prefix_length__(x,p+1)
        else:
            return p
    
    def get_prefix_length(self,x):
        return self.__get_prefix_length__(x,0)
    
    def get_prefix_levels(self,x):
        prefix_len=self.get_prefix_length(x)
        result=None
        
        if prefix_len>=2:
            old_levels=self.get(x[0:prefix_len])
            result=old_levels[0:prefix_len]
        
        return result
    
    def __delete_subtrees__(self,maxlen):
        child_node=self
        while child_node.value is None:                
            child_node=next(iter(child_node.children.values()))
        self.value=child_node.value[0:maxlen]
        self.children={}
    
    def __clean_subtrees__(self,maxlen,p):
        if p<maxlen:
            for child in self.children.values():
                child.__clean_subtrees__(maxlen,p+1)
            
        else:
            self.__delete_subtrees__(maxlen)
            
    def clean_subtrees(self,maxlen):
        self.__clean_subtrees__(maxlen, 0)
        
class WorldState():
    '''
    Represents a state of the world at some time point
    '''     
    def __init__(self):
        self.instances=[]
        self.icpi=None
        self.instance_counts={}
        self.prevalent_pairs=set()
        self.local_candidates=set()
        self.saved_trees=TrieNode()#SortedDict()
        self.local_candidates_divided=[]

    def count_instances(self):
        mydict={}
        for x in self.instances:
            mydict[x.id.feature]=mydict.get(x.id.feature,0)+1
                    
        return mydict
        
    def add_instance(self,instance):
        self.instances.append(instance)
    
    def remove_instance(self,myid):
        for x in range(0,len(self.instances)):
            if self.instances[x].id==myid:
                del self.instances[x]
    
    def build_local_ICPI(self,maxdist):
        
        neighbors=planesweep(self, maxdist)
        
        self.icpi=ICPITree()
        self.icpi.build(neighbors)
    
    def find_prevalent_pairs(self,minprev):
        self.instance_counts=self.count_instances()
        self.prevalent_pairs=self.icpi.prevalent_pairs(minprev, self.instance_counts)
    
    def build_local_candidates(self):
        g=Graph(self.prevalent_pairs)                
        self.local_candidates=g.bron_kerbosch()
    
        lenc=max(len(c) for c in self.local_candidates)
        self.local_candidates_divided=[None for _ in range(lenc+1)]
        while lenc>=2:
            self.local_candidates_divided[lenc]=set(c for c in self.local_candidates if len(c)==lenc)
            lenc=lenc-1
            
            
    
    def remove_prevalent_pairs(self,forbiddenpairs):
        self.prevalent_pairs.difference_update(forbiddenpairs)                
    
    def clean_memory(self,maxlen):
        self.saved_trees.clean_subtrees(maxlen)
    
    def compute_prevalence(self,candidate,save_trees=False):
        ordered_candidate=sorted(list(candidate))
        
        levels=None        
        
        already_built=0
        
        if save_trees:
            tup=tuple(ordered_candidate)
            #keys=self.saved_trees.keys();
            #i = bisect_left(keys, tup)
            #if i:
            #    levels=self.saved_trees[keys[i]]
            #for x in self.saved_trees.keys():                
            #    if x[0:len(ordered_candidate)]==tup:
            #        levels=self.saved_trees[x]                    
            #        break
            
            levels=self.saved_trees.get(tup)
            if levels is None:                
                levels=self.saved_trees.get_prefix_levels(tup)
                if levels is not None:
                    already_built=len(levels)                    
                    levels.extend([] for _ in range(len(ordered_candidate)-already_built))
            else:
                already_built=len(levels) 
                    
                
        
        if levels is None:
            #initialize instance tree  
            levels=[[] for _ in range(len(ordered_candidate))] #make list of lists
        
            #generate first two levels of the instance tree
            p1=ordered_candidate[0]
            p2=ordered_candidate[1] 
            for k in filterfalse(lambda x:not (x[0].feature==p1 and x[1]==p2),self.icpi.dictionary):
                levels[0].append((k[0],None,False))
                pointer=len(levels[0])-1
                for x in self.icpi.find_neighbors(k[0], k[1]):
                    levels[1].append((x,pointer,False))
            
            already_built=2
        
        if already_built<len(ordered_candidate):
            #use icpi tree to build subsequent levels
            
            empty_suffix_len=0;
            for pos in range(already_built,len(ordered_candidate)):
                for pose in range(0,len(levels[pos-1])):
                    e=levels[pos-1][pose]
                    neighbor_set=set(self.icpi.find_neighbors(e[0], ordered_candidate[pos]))
                    parent=e[1]
                    for prev_pos in range(pos-2,-1,-1):
                        fid=levels[prev_pos][parent][0]
                        parent=levels[prev_pos][parent][1]
                        neighbor_set.intersection_update(self.icpi.find_neighbors(fid, ordered_candidate[pos]))
                                    
                    levels[pos].extend((x,pose,False) for x in neighbor_set)
                if len(levels[pos])==0:
                    empty_suffix_len=len(ordered_candidate)-pos
                    break

            #empty_suffice_len=0
            #for x in range(len(levels)-1,-1,-1):
            #    if len(levels[x])==0:
            #        empty_suffix_len=empty_suffix_len+1
            #    else:
            #        break
                    
            if save_trees:                    
                self.saved_trees.insert(tuple(ordered_candidate[0:len(ordered_candidate)-empty_suffix_len]),levels[0:len(ordered_candidate)-empty_suffix_len])
                #self.saved_trees[tuple(ordered_candidate)]=levels
                
            if empty_suffix_len>0:
                return 0
                
        #compute candidate prevalence by finding unique feature instances included in collocation candidate
        unique_fids=[set() for _ in range(len(ordered_candidate))]
        
        for x in levels[len(ordered_candidate)-1]:
            unique_fids[len(ordered_candidate)-1].add(x[0])
            parent=x[1]
            for level in range(len(ordered_candidate)-2,-1,-1):
                (fid,new_parent,visited)=levels[level][parent]
                if (visited):
                    break                
                levels[level][parent]=(fid,new_parent,True)                
                parent=new_parent
                unique_fids[level].add(fid)
        
        if save_trees:
            for x in range(0,len(ordered_candidate)):
                for y in range(0,len(levels[x])):
                    levels[x][y]=(levels[x][y][0],levels[x][y][1],False)
            #for x in levels[len(ordered_candidate)-1]:                
            #    parent=x[1]
            #    for level in range(len(ordered_candidate)-2,-1,-1):
            #        (fid,new_parent,visited)=levels[level][parent]
            #        if not visited:
            #            break                
            #        levels[level][parent]=(fid,new_parent,False)                
            #        parent=new_parent
                    
        
        return min(len(unique_fids[p])/self.instance_counts.get(ordered_candidate[p],1) for p in range(len(ordered_candidate)))

class InputDataset():
    '''
    Represents multiple states of the world in multiple time moments
    '''
    
    def __init__(self):
        self.states={} 
    
    def __row2Instance__(self,row):
        return FeatureInstance(Position(float(row[3]),float(row[4])),FeatureInstanceIdentifier(int(row[1]),int(row[2])))
    
    def loadFromFile(self,filename):
        
        self.states={}
        
        with open(filename,newline='') as csvfile:
            csvreader = csv.reader(csvfile)
            for row in csvreader:
                time=int(row[0])
                if not time in self.states:
                    self.states[time]=WorldState()
                self.states[time].add_instance(self.__row2Instance__(row))
                