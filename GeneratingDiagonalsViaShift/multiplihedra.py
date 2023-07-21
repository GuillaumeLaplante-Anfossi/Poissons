'''
Code for comparing SU and LA diagonal of multiplehdra.
The multipledra filter is code of Guillaume's from,
https://cocalc.com/projects/4ffaeb4c-c563-4526-8bbc-3d75298931c9/files/Multiplihedron_diagonal.ipynb#id=929148
'''
from itertools import combinations

from diagonals_via_shift import SU_diag, LA_diag


def ordered_partitions(a, k):
    '''Splits a into k ordered partitions'''
    n = len(a)
    assert 1 <= k <= n, (n, k)

    def split_at(js):
        i = 0

        for j in js:
            yield a[i:j]
            i = j

        yield a[i:]

    for separations in combinations(range(1, n), k - 1):
        yield list(split_at(separations))

def all_to_right(element,j):
    m = []
    for x in element[j+1:]:
        m = m + x
    return m

def op_mult_filter(n, element):
    '''
    return: True if the ordered partition hits, False if not
    '''
    for j in range(len(element)-1):#for each block in the ordered partition / chosen face
        if len(element[j])==1 or (n in element[j]):
            continue

        que = all_to_right(element,j) 

        OSP=list(ordered_partitions(element[j],2)) #all ordered partitions of this block into two pieces
        for k in range(len(OSP)):#for each block bipartition
            for x in que: #for each element in the queue
                if min(element[j]) < x and x < max(element[j]):
                #if ((max(OSP[k][0]) < x) and (x < min(OSP[k][1]))) or ((max(OSP[k][1]) < x) and (x < min(OSP[k][0]))):
                    #print(OSP[k][0],OSP[k][1])
                    return False

    return True


def LA_Mult_diag(n):
    valid = []
    for p in LA_diag(n):
        if op_mult_filter(n, p[0]) and op_mult_filter(n, p[1]):
            valid.append(p)

    return valid

def SU_Mult_diag(n):
    valid = []
    for p in SU_diag(n):
        if op_mult_filter(n, p[0]) and op_mult_filter(n, p[1]):
            valid.append(p)

    return valid

def compare(n):
    LA = sorted(LA_Mult_diag(n))
    SU = sorted(SU_Mult_diag(n))

    LA_only = []
    SU_only = []
    shared = []
    for x in LA:
        if x in SU:
            shared.append(x)
        else:
            LA_only.append(x)
    for x in SU:    
        if x not in LA:
            SU_only.append(x)

    return LA, LA_only, shared, SU_only, SU


def iso1(p):
    '''
    Applies the iso t(s x s)
    '''
    return [p[1][::-1],p[0][::-1]]

def iso2(p,n):
    '''
    Applies the iso t(r x r)
    '''
    def r(op,n):
        rev = []
        for block in op:
            revblock = []
            for x in block:
                revblock.append(n-x+1)
            rev.append(sorted(revblock))
        return rev

    return [r(p[1],n),r(p[0],n)]

if __name__ == '__main__':

    print(iso2([[[1,3],[2,4]],[[1,3],[2,4]]],4))

    for n in range(1,7):
        print("\nn={}".format(n))
        LA, LA_only, shared, SU_only, SU = compare(n)

        print("|LA|={}, |LA only|={}, |shared|={}, |SU only|={}, |SU|={}".format(
            len(LA), len(LA_only), len(shared), len(SU_only), len(SU)))

        '''
        #ILA = [iso1(x) for x in LA]
        ILA = [iso2(x,n) for x in LA]
        
        ILA_only = []
        SU_only = []
        shared = []
        for x in ILA:
            if x in SU:
                shared.append(x)
            else:
                ILA_only.append(x)
        for x in SU:
            if x not in ILA:
                SU_only.append(x)

        print(len(ILA_only),len(shared),len(SU_only))
        '''