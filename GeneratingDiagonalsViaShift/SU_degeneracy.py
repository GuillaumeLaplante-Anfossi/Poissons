'''
SU degeneracy projections, 
using Proposition 6 of DIAGONALS ON THE PERMUTAHEDRA, MULTIPLIHEDRA AND ASSOCIAHEDRA.
'''
from itertools import combinations

from diagonals_via_shift import SU_diag, LA_diag


def exceptional(A,k):
    '''
    return: True if the subset A_k is exceptional, False if not
    '''
    minAk = min(A[k])
    maxAk = max(A[k])

    for x in A[k+1:]:
        for aij in x:
            if minAk < aij and aij < maxAk:
                return True

    return False

def Mult_filter(A):
    for j in range(len(A)-1):
        minAj = min(A[j])
        if (min(A[j]) > min([min(x) for x in A[j+1:]])) and exceptional(A,j):
            return False

    return True

def LA_Mult_diag(n):
    valid = []
    for p in LA_diag(n):
        if Mult_filter(p[0]) and Mult_filter(p[1]):
            valid.append(p)

    return valid

def SU_Mult_diag(n):
    valid = []
    for p in SU_diag(n):
        if Mult_filter(p[0]) and Mult_filter(p[1]):
            valid.append(p)

    return valid

def Assoc_filter(A):
    for k in range(len(A)-1):
        if exceptional(A,k):
            return False

    return True

def LA_Assoc_diag(n):
    valid = []
    for p in LA_diag(n):
        if Assoc_filter(p[0]) and Assoc_filter(p[1]):
            valid.append(p)

    return valid

def SU_Assoc_diag(n):
    valid = []
    for p in SU_diag(n):
        if Assoc_filter(p[0]) and Assoc_filter(p[1]):
            valid.append(p)

    return valid

def compare_Mult(n):
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

def compare_Assoc(n):
    LA = sorted(LA_Assoc_diag(n))
    SU = sorted(SU_Assoc_diag(n))

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

def iso1(p):
    '''
    Applies the iso t(s x s)
    '''
    return [p[1][::-1],p[0][::-1]]

if __name__ == '__main__':

    print("\nAssoc")
    for n in range(1,6):
        print("\nn={}".format(n))
        LA, LA_only, shared, SU_only, SU = compare_Assoc(n)

        print("|LA|={}, |LA only|={}, |shared|={}, |SU only|={}, |SU|={}".format(
            len(LA), len(LA_only), len(shared), len(SU_only), len(SU)))

    print("\nMult")
    for n in range(1,8):
        print("\nn={}".format(n))
        LA, LA_only, shared, SU_only, SU = compare_Mult(n)

        print("|LA|={}, |LA only|={}, |shared|={}, |SU only|={}, |SU|={}".format(
            len(LA), len(LA_only), len(shared), len(SU_only), len(SU)))

