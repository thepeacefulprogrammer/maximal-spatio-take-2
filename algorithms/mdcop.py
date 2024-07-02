'''
Created on 22.12.2016

@author: witek
'''

from itertools import takewhile
from itertools import groupby
from bisect import  bisect_left

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    #raise ValueError
    return -1

def join(f,verbose=0):
    c=[] #result list with candidates
    p=[] #result list with positions of candidate "fathers"
    
    if len(f)==0:
        return c,p
    
    prefixlen=len(f[0])-1
    for pos in range(0,len(f)):
        if verbose>1:
            print("a: ",f[pos])
        prefix=f[pos][0:prefixlen]
        if verbose>1:
            print("prefix: ",prefix)
        for mate in takewhile(lambda x:prefix==x[0:prefixlen],(f[i] for i in range(pos+1,len(f)))):
            if verbose>1:
                print("mate: ", mate)
                
            candidate=f[pos]+(mate[prefixlen],)  
            
            success=True
            for sc in range(0,len(candidate)-2):
                subcandidate=tuple(candidate[i] for i in range(0,len(candidate)) if i!=sc)
                if verbose>1:
                    print("subcandidate: ",subcandidate)
                    
                if (index(f,subcandidate)==-1):
                    success=False
                    break
                
            if success:
                c.append(candidate)
                p.append(pos)
            
            
            
    return (c,p)



def upper_bound(f,c,total):    
    return (total-(c-f))/total

def upper_bound_frames(f,c,total):    
    return (total-(c-f))

def build_instances(extending_feature,prefix_instances,icpi):
    result=[]
    for instance in prefix_instances:
        common_neighbors=set(icpi.find_neighbors(instance[0],extending_feature))
        for p in range(1,len(instance)):
            common_neighbors.intersection_update(icpi.find_neighbors(instance[p],extending_feature))
            if len(common_neighbors)==0:
                break
            
        for x in sorted(common_neighbors):            
            result.append(instance+(x,))
    
    return result

def mdcop(input,maxdist,minprev,minfreq,verbose=0):
    
    candidates={}
    colocations={}
    
    minfreq_frames=minfreq*len(input.states)
    
    
    #Find prevalent pairs in each time moment, remove those that can't be frequent on the fly
    forbidden_pairs=set()
    frequencies={}
    count=0 
    
    for statenum in input.states.keys():
        count=count+1
        state=input.states[statenum]
        state.instance_counts=state.count_instances()
        state.build_local_ICPI(maxdist)
        (pairs,instances)=state.icpi.prevalent_pairs_with_instances(minprev,state.instance_counts,forbidden_pairs)
        colocations[(2,statenum)]=(pairs,instances)
        for p in pairs:
            frequencies[p]=frequencies.get(p,0)+1
            if upper_bound_frames(frequencies[p], count, len(input.states))<minfreq_frames:
                forbidden_pairs.add(p)
                del frequencies[p]
    
    #Compute final frequencies
    for p in frequencies.keys():            
        if upper_bound_frames(frequencies[p], count, len(input.states))<minfreq_frames:
            forbidden_pairs.add(p)
                
    #remove non frequent pairs
    for statenum in input.states.keys():
        (pairs,instances)=colocations[(2,statenum)]        
        new_pairs=[]
        new_instances=[]
        for i in range(0,len(pairs)):
            if not (pairs[i] in forbidden_pairs):
                new_pairs.append(pairs[i])
                new_instances.append(instances[i])
                
        colocations[(2,statenum)]=(new_pairs,new_instances)
    
    #Build longer co-locations
    clen=2
    while True:
        clen=clen+1
        
        #Build candidates by joining shorter co-locations
        newcandidates=0;
        for statenum in input.states.keys():
            candidates[(clen,statenum)]=join(colocations[(clen-1,statenum)][0])
            newcandidates=newcandidates+len(candidates[(clen,statenum)][0])
        
        if newcandidates==0:
            break
                
        #Verify candidates
        count=0
        frequencies={}
        forbidden_candidates=set()

        for statenum in input.states.keys():
            count=count+1
            (state_candidates,positions)=candidates[(clen,statenum)]
            local_colocations=[]
            local_instances=[]
            for cnum in range(0,len(state_candidates)):
                if state_candidates[cnum] in forbidden_candidates:
                    continue
                
                
                
                instances=build_instances(state_candidates[cnum][len(state_candidates[cnum])-1],colocations[(clen-1,statenum)][1][positions[cnum]],input.states[statenum].icpi)                
                
                
                #W jednej linijce liczę całe prevalence. Jestem z siebie dumny :D
                prevalence=min((len({instances[i][x] for i in range(0,len(instances))})/input.states[statenum].instance_counts[state_candidates[cnum][x]] for x in range(0,clen)))
                
                                
                if prevalence>=minprev:
                    frequencies[state_candidates[cnum]]=frequencies.get(state_candidates[cnum],0)+1
                    local_colocations.append(state_candidates[cnum])
                    local_instances.append(instances)
                    
                                
                if upper_bound_frames(frequencies.get(state_candidates[cnum],0), count, len(input.states))<minfreq_frames:   
                    forbidden_candidates.add(state_candidates[cnum])
                    frequencies.pop(state_candidates[cnum],None) #Use pop to avoid keyerror
            
            colocations[(clen,statenum)]=(local_colocations,local_instances)
        
        #Compute final frequencies
        for c in frequencies.keys():
            if upper_bound_frames(frequencies[c], count, len(input.states))<minfreq_frames:   
                    forbidden_candidates.add(c)
        
        # remove non_frequent co-locations
        new_colocation_count=0;
        for statenum in input.states.keys():
            (cols,instances)=colocations[(clen,statenum)]        
            new_cols=[]
            new_instances=[]
            for i in range(0,len(cols)):
                if not (cols[i] in forbidden_candidates):
                    new_cols.append(cols[i])
                    new_instances.append(instances[i])
            
            new_colocation_count=new_colocation_count+len(new_cols)
            colocations[(clen,statenum)]=(new_cols,new_instances)
            
        if new_colocation_count==0:
            break;
        
    final_colocations=set()
    for col in colocations.values():
        final_colocations.update(col[0])
    
    return final_colocations

def is_subtuple(a,b):
    if len(a)==0:
        return True
    if len(b)==0:
        return False
    pa=0
    pb=0
    while pa<len(a) and pb<len(b):
        while pb<len(b) and a[pa]>b[pb]: #WARNING! partial evaluation of boolean equations is required for this to work correctly
            pb=pb+1
        if pb==len(b):
            return False            
        if a[pa]!=b[pb]:
            return False
        pa=pa+1  
    return True

def find_maximal(colocs):
    sorted_colocs=sorted(colocs,key=lambda x:len(x),reverse=True)
    results=[]
    for k,g in groupby(sorted_colocs,lambda x:len(x)):
        temp=[]
        for c in g:
            if not any(is_subtuple(c, x) for x in results):
                temp.append(c)
        results.extend(temp)    
    return results