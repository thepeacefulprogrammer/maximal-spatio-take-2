'''
Created on 30.11.2016

@author: witek
'''
from functools import reduce
from time import perf_counter

def output_prevalent_pairs(input,all_pairs,verbose):
    if (verbose>1):
        print("Prevalent pairs for each time moment:")
        for time in input.states.keys():
            print("Time t=",time)
            state=input.states[time]
            print(state.prevalent_pairs)
        print("All pairs:")
        print(all_pairs)

def output_nonfrequent_pairs(forbidden_pairs,verbose):
    if (verbose>1):
        print("Non frequent pairs")
        print(forbidden_pairs)

def output_size2_MDCOPS(input,verbose):
    if (verbose>1):
        print("Prevalent and frequent pairs for each time moment:")
        for time in input.states.keys():
            print("Time t=",time)
            state=input.states[time]
            print(state.prevalent_pairs)

def output_global_candidates(global_candidates,verbose):
    if (verbose>1):
        print("Pruned global candidates:")
        print(global_candidates)

def output_local_and_global_candidate_sets(input,initial_global_candidates,verbose):
    if (verbose>1):
        print("Local candidates:")
        for time in input.states.keys():
            print("Time t=",time)
            state=input.states[time]
            print(state.local_candidates)
            
        print("Initial global candidates (including non frequent ones)")
        print(initial_global_candidates)
    

def output_results_length_cache(maxlen,result,cache,verbose):
    if (verbose>1):
        print("Current length ",maxlen)
        print("\tCurrent result set: ", result)
        print("\tCurrent cache: ", cache)    

def output_step_time(step,time,verbose):    
    if (verbose>0):    
        print("Step",step,"time",time)

def output_iteration_time(maxlen,time,verbose):
    if verbose>0:
        print("Main, length ",maxlen,":",time)
    
def output_cache_stats(cache,cachetime, cachehit,cachemiss,verbose):
    if (verbose>0):
        print("Cache time", cachetime, "cache hit ",cachehit, "cachemiss ", cachemiss, "cachesize", reduce(lambda x,y:x+y,(len(x) for x in cache.values())))
    
def output_glocal_candidate_stats(maxlen,global_candidates,verbose):
    if (verbose>0):
        if maxlen-1 in global_candidates:
            print("Candidate size ",maxlen-1," count ",len(global_candidates[maxlen-1]))
        else:
            print("Candidate size ", maxlen - 1, " count 0")
    
def output_iteration_state(maxlen,result,cache,new_candidates,verbose):
    if verbose > 1:
        print("Current length ",maxlen)
        print("\tCurrent result set: ", result)
        print("\tCurrent cache: ", cache)
        print("\tNew candidates added:", new_candidates)
        
#Finds prevalent pairs in each time moment (stored as a field prevalent_pairs in the WorldState objects)
#as well as a set of all prevalent pairs (returned as a function value)
def find_prevalent_pairs(input,maxdist,minprev):
    #Declare a set which will contain all frequent pairs
    all_pairs=set()
    
    #Find prevalent pairs in each time moment and put them all into one set
    for state in input.states.values():
        #Build an iCPI tree for the time moment represented by a state
        #iCPI tree is built using plane sweep
        state.build_local_ICPI(maxdist)        
        #Convert neighborhood information from the iCPI tree into a set of prevalent pairs (local variable in the state object) 
        state.find_prevalent_pairs(minprev)
        #Compute set sum of the up-to-date found pairs with the pairs found for the current state
        #Note that duplicates are removed since all_pairs is a set
        all_pairs.update(state.prevalent_pairs)
    
    #Return the found set of frequent pairs
    return all_pairs

#Finds pairs that are prevalent in some time moments, but are infrequent (returned as a function value)
def find_nonfrequent_pairs(input,all_pairs,minfreq_it):
    #Declare a set for storing infrequent pairs
    forbidden_pairs=set()
    
    #For all pair prevalent at some time moments 
    for pair in all_pairs:
        #compute frequency of a pair expressed as a number of time moments in which it is prevalent
        #freq=sum(1 for x in input.states.values() if pair in x.prevalent_pairs)/len(input.states)
        freq_it=sum(1 for x in input.states.values() if pair in x.prevalent_pairs)
        #If the computed frequency is smaller than the minfreq threshold then add the pair to the result set
        if freq_it<minfreq_it:
            forbidden_pairs.add(pair)
    
    #Return the found set of infrequent pairs
    return forbidden_pairs


#Modifies states of each time moment by removing nonfrequent pairs from prevalent_pairs sets
def remove_nonfrequent_pairs(input,forbidden_pairs):
    #Update prevalent_pairs in all time moments to store only time prevalent pairs
    for state in input.states.values():
        state.remove_prevalent_pairs(forbidden_pairs)        


#Generates local candidate sets and their union: a global candidate set
def build_local_and_global_candidate_sets(input):
    #Declare a map for storing global candidates
    #Each key will represent size of a candidate set. The corresponding value will be a set of candidate sets of proper size 
    initial_global_candidates={}
    
    #For each time moment...
    for state in input.states.values():
        #...build a set of local candidates based of frequent prevalent pairs
        state.build_local_candidates()
        #For each of the local candidates insert it into an appropriate set in the initial_global_candidates map 
        for c in state.local_candidates:
            #If the appropriate entry in the initial_global_candidates map is not created yet then create it/ 
            if not len(c) in initial_global_candidates: 
                initial_global_candidates[len(c)]=set()
            #Insert the candidate c into appropriate set of candidates
            initial_global_candidates[len(c)].update({c})
    
    #Return the set of global_candidates
    return initial_global_candidates


   
#Modifies initial global candidate set by removing candidates that are infrequent.
#This is done by checking whether a candidate is a subset of enough local candidates to be frequent
#In case a candidate is not frequent among other candidates, it is replaced by all of its size-1 subsets   
def prune_candidate_set(input,initial_global_candidates,minfreq_it):
    #Declare map for storing pruned candidate sets
    #Each key will represent size of a candidate set. The corresponding value will be a set of candidate sets of proper size 
    global_candidates={}
    
    #As long as there are still sets of candidates in the initial global candidates map...
    while len(initial_global_candidates)>0:
        #Get maximal candidate length yet unprocessed
        maxlen=max(initial_global_candidates.keys())
        #Initialize an entry in the global_candidates map for candidates of size maxlen
        global_candidates[maxlen]=set()
        
        #Extract and remove largest candidates from the initial global candidates map
        candidates=initial_global_candidates[maxlen]        
        initial_global_candidates.pop(maxlen)
        
        #For each of the extracted candidates
        for candidate in candidates:
            
            #Check if superset is already added, if so, omit the candidate
            remove=False
            for k in global_candidates.keys():
                if any(candidate<=a for a in global_candidates[k] if k>maxlen):
                    remove=True
                    break                     
                    
             
            if (remove):
                continue
            
            #Compute frequency of a candidate
            #freq=reduce(lambda x,y:x+y,(1 for x in input.states.values() if any(candidate<=a for a in x.local_candidates)) )/len(input.states)
            freq_it=reduce(lambda x,y:x+y,(1 for x in input.states.values() if any(candidate<=a for a in x.local_candidates)) )
            
            #if candidate is potentially frequent then add it to filtered global candidates
            #otherwise test all of its subsets by adding them to the initial global candidates         
            if freq_it>=minfreq_it:
                global_candidates[maxlen].add(candidate)
            else:
                if len(candidate)>2: # Do not add subsets of size 2 sets
                    for item in candidate:
                        if not maxlen-1 in initial_global_candidates:
                            initial_global_candidates[maxlen-1]=set()
                        initial_global_candidates[maxlen-1].add(candidate.difference({item}));
    
        #If no global candidates of length maxlen are potentially frequent then remove the corresponding entry from the global candidates map
        if len(global_candidates[maxlen])==0:
            global_candidates.pop(maxlen,0)
        
        #If size 2 candidates are have been processed, than we can end this loop
        if maxlen==2:
            break
        
        
    #Return the pruned global candidates map
    return global_candidates


#Procedure iterates over every instance tree in every time moment and deletes every tree level at maxlen depth.
def clean_tree_memory(input,maxlen):
    for statenum in input.states: # for each of the time moments
        state=input.states[statenum] 
        state.clean_memory(maxlen) #delete every tree level at maxlen depth

#Procedure finds frequent and prevalent pairs.
#This procedure is used instead of normal mining algorithm, since size 2 prevalent and frequent pairs are found at the start of the algorithm
def process_size2_candidates(all_pairs,forbidden_pairs,global_candidates):
    #Remove non frequent pairs from the set of all prevalent pairs
    #This could have been done earlier, but was not needed and might not be needed at all
    all_pairs.difference_update(forbidden_pairs)
    #Convert tuple representation to frozen set representation of the pairs
    #TODO: Think on how to avoid converting form 
    all_pairs={frozenset({x,y}) for (x,y) in all_pairs} 
    #Find intersection of the set of frequent and prevalent pairs with the set of global candidates of size 2
    #This yields the set of maximal size 2 MDCOPS
    all_pairs.intersection_update(global_candidates[2])


#Procedure computes frequency of a candidate. Frequency is computed in a form of a set of time moment identifiers in which the candidate is prevalent.
#Consequently, real frequency is in fact the size of the set. The set from_cache stores the time moments in which the prevalence of the candidate was determined by the cache lookup instead of prevalence computation. 
#Also, procedure measures time spent on cache lookups and the number of cache hits and misses. The measurements are returned as a tuple
#
#The arguments are as follows:
#input - the input dataset with all the time moments
#candidate - the candidate for which the frequency is computed
#cache - the cache map
#minfreq_it - minimal frequency threshold - ignored
#minprev threshold
#frequency - result set
#from_cache - result set
#save_trees - true if instance trees should be saved for use in subsequent prevalence computations
def find_frequency_without_estimations(input,candidate,cache,minfreq_it,minprev,frequency,from_cache,save_trees):
    cachetime=0;
    cachehit=0;
    cachemiss=0;
    
    for statenum in input.states: # for each of the time moments
        state=input.states[statenum]
                                                                        
        #check in cache if the candidate is a subset of a locally  prevalent candidate
        t3=perf_counter()
        incache=any(candidate<=a for a in cache[statenum])
        t4=perf_counter()
        cachetime=cachetime+t4-t3
                                
        if incache:
            cachehit=cachehit+1
            frequency.append(statenum) #remember the time moment at which the candidate is prevalent
            from_cache.add(statenum) #remember the time moment at which the cache allowed to determine that the candidate is prevalent 
        else:
            cachemiss=cachemiss+1
            #if not cached, then compute candidate's prevalence in the current time moment                
            prevalence=state.compute_prevalence(candidate,save_trees);
            if prevalence>=minprev:
                frequency.append(statenum)  #remember the time moment at which the candidate is  prevalent
    
    return (cachetime,cachehit,cachemiss)

#Procedure computes frequency of a candidate. Frequency is computed in a form of a set of time moment identifiers in which the candidate is prevalent.
#To reduce the time of computations, the procedure computes upper and lower bounds on the frequency based on the known results and stops computations 
#if the upper bound is smaller than minfreq threshold or lower bound is higher than minfreq threshold.
#Consequently, lower bound of frequency is in fact the size of the set. The set from_cache stores the time moments in which the prevalence of the 
#candidate was determined by the cache lookup instead of prevalence computation. 
#Also, procedure measures time spent on cache lookups and the number of cache hits and misses. The measurements are returned as a tuple
#
#The arguments are as follows:
#input - the input dataset with all the time moments
#candidate - the candidate for which the frequency is computed
#cache - the cache map
#minfreq_it - minimal frequency threshold
#minprev threshold
#frequency - result set
#from_cache - result set
#save_trees - true if instance trees should be saved for use in subsequent prevalence computations
def find_frequency_with_estimations(input,candidate,cache,minfreq_it,minprev,frequency,from_cache,save_trees):
    cachetime=0;
    cachehit=0;
    cachemiss=0;
    
    #Initialize upper and lower bounds
    #frequ=1
    #freql=0
    frequ_it=len(input.states)
    freql_it=0
    
    #The first loop is based on the cache only.
    #It iterates over every time moment and updates the lower bound based solely on the information from caches.
    #Upper bound cannot be computed since it is based on the information in which time moments candidate IS NOT prevalent
    #whereas cache can only give the information that the candidate IS prevalent.
    for statenum in input.states: # for each of the time moments
                                                                        
        #check in cache if the candidate is a subset of a locally prevalent candidate
        t3=perf_counter()
        incache=any(candidate<=a for a in cache[statenum])
        t4=perf_counter()
        cachetime=cachetime+t4-t3
        
        
        if incache:
            cachehit=cachehit+1
            frequency.append(statenum) #remember the time moment at which the candidate is prevalent                        
            from_cache.add(statenum) #remember the time moment at which the cache allowed to determine that the candidate is prevalent
            #freql=len(frequency)/len(input.states)
            freql_it=len(frequency) #update lower bound
                                       
        #if the lower bound is greater than minfreq threshold no further computations are necesary
        if freql_it>=minfreq_it:
            break
        
    #if the lower bound is smaller than minfreq threshold we need real prevalence computations
    if freql_it<minfreq_it:
        j=len(frequency) #The number of time moments the candidate has been checked up to now. 
        
        time_moments=set(input.states.keys()).difference(set(frequency)) #The set of time_moments in which cache did not give answer
        
        for statenum in time_moments:
            state=input.states[statenum]
            j=j+1
            
            #prefilter based on local candidates - a candidate must be a subset of a local candidate in order to be prevalent at a given time moment
            #incandidates=any(candidate<=a for a in input.states[statenum].local_candidates)
            
            incandidates=None
            lenc=len(input.states[statenum].local_candidates_divided)-1
            while lenc>=len(candidate):
                incandidates=any(candidate<=a for a in input.states[statenum].local_candidates_divided[lenc])
                if incandidates:
                    break
                lenc=lenc-1
            
            if incandidates:                   
                cachemiss=cachemiss+1
                #Compute candidate's prevalence in the current time moment                
                prevalence=state.compute_prevalence(candidate,save_trees)
                
                #If the candidate is prevalent in the current time moment
                if prevalence>=minprev:
                    frequency.append(statenum)  #remember the time moment at which the candidate is spatially prevalent
                    #Update lower bound
                    #freql=len(frequency)/len(input.states)
                    freql_it=len(frequency)
                    
                '''
                else:
                    #prune local candidates
                    new_local_candidates=set();
                    candidate_for_removal=set();
                    
                    #brzydkie i tymczasowe
                    if type(input.states[statenum].local_candidates) is not set:
                        input.states[statenum].local_candidates=set(input.states[statenum].local_candidates)
                        
                    for c in input.states[statenum].local_candidates:
                        if candidate<=c:
                            candidate_for_removal.update({c})                        
                            for x in candidate:
                                new_local_candidates.update({c.difference({x})})
                    
                    input.states[statenum].local_candidates.difference_update(candidate_for_removal)
                    input.states[statenum].local_candidates.update(new_local_candidates)
                            
                    #input.states[statenum].local_candidates=new_local_candidates
                            
                    #input.states[statenum].local_candidates.difference_update(candidate_for_removal);
                    #input.states[statenum].local_candidates.update(new_local_candidates);
                '''
                    
            #Update upper bound
            #frequ=1-((j-len(frequency))/len(input.states))
            frequ_it=len(input.states)-(j-len(frequency))
            
            #Break the loop if no more computations are necessary
            if frequ_it<minfreq_it or freql_it>=minfreq_it:
                break
            
    return (cachetime,cachehit,cachemiss)


#Add new candidates to the global_candidates map. Only new candidates that are not subsets of found results are added
def update_global_candidates(global_candidates,new_candidates,maxlen,result):
    #check if new candidates are not subsets of current results
    for candidate in new_candidates:            
        if not any(candidate<=a for a in result): #If new candidate is not a subset of any result then...
            if not maxlen-1 in global_candidates:
                global_candidates[maxlen-1]=set()
            global_candidates[maxlen-1].add(candidate) #... add it to global_candidates
          
def maxspatiotempcolloc(input,maxdist,minprev,minfreq,predictions=False,save_trees=0,verbose=0,clean_trees=0):
    
    result=[]
    
    minfreq_it=minfreq*len(input.states)
        
    start=perf_counter()
        
    #Find prevalent pairs in each time moment as well as a set of distinct pairs that are prevalent at some time moment
    all_pairs=find_prevalent_pairs(input, maxdist, minprev)     
    output_prevalent_pairs(input, all_pairs,verbose)
                            
    #Find and remove non-time prevalent pairs 
    forbidden_pairs=find_nonfrequent_pairs(input, all_pairs,minfreq_it)    
    output_nonfrequent_pairs(forbidden_pairs,verbose)        
    remove_nonfrequent_pairs(input,forbidden_pairs)
    output_size2_MDCOPS(input,verbose)
    
    end=perf_counter()
    
    output_step_time(1,end-start,verbose)
         
    start=perf_counter()
        
    #Find initial local and global candidates by finding maximal cliques in graphs built from prevalent pairs in each time moment (local candidates)
    #and then finding a set of distinct local candidates
    initial_global_candidates=build_local_and_global_candidate_sets(input)
    
    end=perf_counter()
    
    output_step_time(2,end-start,verbose)
    output_local_and_global_candidate_sets(input, initial_global_candidates,verbose)
    
    start=perf_counter()
    
    #remove candidates that cannot be time-prevalent    
    global_candidates=prune_candidate_set(input,initial_global_candidates,minfreq_it)
    end=perf_counter()
    
    output_step_time(3,end-start,verbose)
    output_global_candidates(global_candidates,verbose)
    
        
    #Check which candidates are prevalent - i.e. mine maximal MDCOPS        
    
    #Initialize cache 
    cache={x:set() for x in input.states.keys()}
            
    #While there are some candidates left
    while len(global_candidates)>0:
        
        start=perf_counter()
        
        #Find the size of the largest candidates
        #maxlen=max(len(x) for x in global_candidates) # get max candidate length
        maxlen=max(global_candidates.keys())
        
        #Perform optional tree cleaning if the tree saving options is turned on        
        if clean_trees and save_trees:
            clean_tree_memory(input,maxlen)
                
        #For candidates of size 2, frequent and prevalent pairs have already been found. 
        if maxlen==2: 
            #Just find intersection of candidates with all_pairs\forbidden_pairs
            process_size2_candidates(all_pairs, forbidden_pairs, global_candidates)
            #Add the found intersection to the result set
            result.extend(all_pairs)
            output_results_length_cache(maxlen, result, cache, verbose)
                        
        #For candidates of size >2 proceed with normal computations
        else:
            #get longest candidates
            #candidates=[x for x in global_candidates if len(x)==maxlen]
            candidates=global_candidates[maxlen]
            
            #remove longest candidates from global_candidates
            #global_candidates.difference_update(candidates) 
            global_candidates.pop(maxlen)
            
            #The set will contain candidates that are generated while mining the current candidates
            new_candidates=set()
        
            cachetime=0
            cachehit=0
            cachemiss=0
            
            # for each candidate verify whether it is frequent or not
            for candidate in candidates: 
                frequency=[] # time moments at which the candidate is prevalent
                from_cache=set() # the set stores information in which time moments prevalence of the candidate stems from cache lookup instead of prevalence computations 
      
                #Depending whether frequency condition should be checked based on lower/upper bounds on frequency or not
                if predictions:
                    (_cachetime,_cachehit,_cachemiss)=find_frequency_with_estimations(input, candidate, cache, minfreq_it, minprev, frequency, from_cache, save_trees)                        
                else:            
                    (_cachetime,_cachehit,_cachemiss)=find_frequency_without_estimations(input, candidate, cache, minfreq_it, minprev, frequency, from_cache, save_trees)
                
                (cachetime,cachehit,cachemiss)=(cachetime+_cachetime,cachehit+_cachehit,cachemiss+_cachemiss)  
                                
                #freq=len(frequency)/len(input.states) #divide the frequency by the number of time moments to compute the time-prevalence of the candidate
                #Compute frequency of the candidate (if no bounds are used) or lower bound (if bounds are used)
                freq_it=len(frequency) 
                                
                if freq_it>=minfreq_it:
                    result.append(candidate) #if the candidate is frequent add it to the result set
                else:
                    #otherwise, test all of its subsets
                    if maxlen>2: #We generate subsets of sets larger than 2, since mining MDCOPS of size 1 is nonsense
                        for item in candidate:
                            new_candidates.add(candidate.difference({item}))
                        
                    #add the candidate to cache for every time moment at which the candidate is spatially prevalent
                    #add only if maxlen>3 since caching mechanism will not be used for candidates of length 2
                    if maxlen>3:                    
                        for statenum in set(frequency).difference(from_cache):
                            cache[statenum].add(candidate)
                            
            output_cache_stats(cache,cachetime, cachehit,cachemiss,verbose)

            #Add new candidates to the global_candidates map. Only new candidates that are not subsets of found results are added
            update_global_candidates(global_candidates,new_candidates,maxlen,result)                        

            output_glocal_candidate_stats(maxlen,global_candidates,verbose)
            output_iteration_state(maxlen,result,cache,new_candidates,verbose)
            

        end=perf_counter()
  
        output_iteration_time(maxlen, end-start, verbose)          
                        
        #If size 2 candidates have been processed, than we can end this loop
        if maxlen==2:
            break;          
    
    return result