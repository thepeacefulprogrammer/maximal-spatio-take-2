
'''
Created on 18.12.2016

@author: witek
'''

from bitstring import BitArray

class HierarchicalBitmap(object):
    '''
    classdocs
    '''


    def __init__(self, nodesize, myset):
        '''
        Constructor
        '''
        self.nodesize=nodesize
        
        
        #compute height
        maxvalue=max(myset)
        
        h=1
        
        while (maxvalue>nodesize):
            maxvalue=maxvalue//nodesize
            h=h+1
        
        
        # build levels
        
        self.levels=[{} for _ in range(0,h)]
        for v in myset:
            x=v
            for l in range(h-1,-1,-1):
                self.levels[l][x//self.nodesize]=self.levels[l].get(x//self.nodesize,0)|(x%self.nodesize)
                x=x//self.nodesize 
        
        
        self.keys=[sorted(x.keys()) for x in self.levels]        
        
        

    def height(self):
        return len(self.levels)

    
    def __extend_height_(self,height):
        self.hbmp=[{0:1} for _ in range(0,height-self.height())]
        self.keys=[[0] for _ in range(0,height-self.height())]
        


    def __is_superset__(self,other):
        '''
        Checks if the object represents a superset of the other object
        '''
        return other.__is_subset__(self)
                

    def __is_subset__(self,other):
        '''
        Checks if the object represents a subset of the other object
        '''
        
        h1=self.height()
        h2=other.height()
        
        if h1<h2:
            for h in range (0,h2-h1):
                if other.keys[h][0]!=0:
                    return False
                else:
                    if (1&other.levels[h][0]!=1):
                        return False
            b1=h2-h1
            b2=0
            
        if h1>h2:
            for h in range(0,h1-h2):
                if self.keys[h][0]!=0 or len(self.keys[h])>1:
                    return False
                else:
                    if (self.levels[h][0]!=1):
                        return False
        
            b1=0
            b2=h1-h2
        
        for h in range(0,min(h1,h2)):
            p2=0
            for p1 in range(0,len(self.keys[b1+h])):
                k1=self.keys[b1+h][p1]
                while other.keys[b2+h][p2]<k1:
                    p2=p2+1
                    if p2>len(other.keys[b2+h][p2])-1:
                        return False 
                    else:
                        if other.keys[b2+h][p2]>k1:
                            return False
            
                
                if (self.levels[b1+h][k1]&other.levels[b2+h][other.keys[b2+h][p2]])!=self.levels[b1+h][k1]:
                    return False
                
        return True
    
    def __is_equal__(self,other):
        '''
        Checks if the object represents a set equal to the other object
        '''
        h1=self.height()
        h2=other.height()
        
        if h1<h2:
            for h in range (0,h2-h1):
                if other.keys[h][0]!=0:
                    return False
                else:
                    if (other.levels[h][0]!=1):
                        return False
            b1=h2-h1
            b2=0
            
        if h1>h2:
            for h in range(0,h1-h2):
                if self.keys[h][0]!=0 or len(self.keys[h])>1:
                    return False
                else:
                    if (self.levels[h][0]!=1):
                        return False
        
            b1=0
            b2=h1-h2
            
        for h in range(0,min(h1,h2)):
            if (len(self.keys[b1+h])!=len(other.keys[b2+h])):
                return False
                        
            for p1 in range(0,len(self.keys[b1+h])):
                k1=self.keys[b1+h][p1]
                if (other.keys[b2+h][p1]!=k1):
                    return False
                
                if self.levels[b1+h][k1]!=other.levels[b2+h][k1]:
                    return False
        
        return True

    def __set_bit__(self,number):
        '''
        Sets a bit (adds an item)
        '''
        
        
    def __clear_bit__(self,number):
        '''
        Clears a bit (clears an item)
        '''

    def __bitwise_and__(self,number):
        '''
        Performs a bitwise-and of the two bitmaps (computes intersection of the sets)
        '''

    def __bitwise_or__(self,number):
        '''
        Performs a bitwise-or of the two bitmaps (computes a sum of the sets)
        '''

    def __bitwise_xor__(self,number):
        '''
        Performs a bitwise-xor of the two bitmaps
        '''

    def __negate__(self):
        '''
        Negates a Bitmap
        '''
    def __intersection_not_empty__(self,other):
        '''
        Returns true if bitwise-and would return a bitmap filled with zeroes
        '''

    
'''
def test_bit(bmp,bitnum):
    return (bmp&(1<<bitnum))!=0

def hbmp_generator_rek(hbmp,parentnode,basenum,level,maxlevel):
    if level==maxlevel:
        pos=0
        for x in range(0,hbmp.nodelen):
            if (test_bit(parentnode)):
                yield (x,hbmp.levels[level][basenum+pos])
                pos=pos+1
    else:
        pos=0
        total=0
        for x in range(0,hbmp.nodelen):
            if (test_bit(parentnode)):
                total+=1<<x
                parentnode=hbmp.levels[level][basenum+pos]
                for p in hbmp_generator_rek(hbmp,parentnode,basenum+pos,level+1):
                    yield (total+p[0],p[1])
                    
                pos=pos+1    


class HierarchialBitmapIterator(object):

    def __init__(self,hbmp):
        self.hbmp=hbmp
        self.levelbitpos=[0 for _ in range(hbmp.height()-1)]
        self.levelpos=[0 for _ in range(hbmp.height()-1)]
        self.levelnodenum=[0 for _ in range(hbmp.height()-1)]
                
        for x in range(0,hbmp.height()):
            self.levelbitpos[x]=__youngest__(self,hbmp.levels[x][0])
            
        
    def __youngest_bit__(self,bmp):
        tmp=bmp
        x=1
        while x<len(bmp):
            tmp=tmp|(tmp<<x)
            x=x*2
        print(tmp.bin)
    
    def __getmask__(self,s,e):
        p1=(-1<<s)
        p2=~(-1<<(e+1))
        return p1&p2

    def __youngest__(self,bmp,nodelength):
        if bmp==0:
            return -1
        tmp=bmp&-bmp
        
        s=0
        e=nodelength-1
        
        while(s<e):
            m=(s+e)//2
            p1=__getmask__(self,s,m)
            
            if tmp&p1!=0:
                e=m
            else:
                s=m+1
        
        return s
    
    def __iter__(self):
        return self
    
    
    
    def __next__(self):
        nodenum=0
        for h in range(self.hbmp.height()-2,-1,-1):
            #Check if node has any additional bits
            levelpos=self.levelpos[h]
            levelbitpos=self.levelbitpos[h]
            if self.hbmp.levels[h][levelpos]&
        
'''   