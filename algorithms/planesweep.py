'''
Created on 28.11.2016

@author: witek
'''


def planesweep(state,distance):
    instances=sorted(state.instances,key=lambda k:k.pos.x)
    
    result=[]
    
    for i in range(1,len(instances)):
        j=i-1
        a=instances[i]
        
        while (a.pos.x-instances[j].pos.x)<=distance:
            b=instances[j]
            
            if (a.pos.distance(b.pos)<=distance):
                if a.id.feature<b.id.feature:
                    result.append((a.id,b.id))
                if a.id.feature>b.id.feature:
                    result.append((b.id,a.id))
                if a.id.feature==b.id.feature:
                    pass
            
            j-=1
            
            if j<0:
                break
            
    return result

