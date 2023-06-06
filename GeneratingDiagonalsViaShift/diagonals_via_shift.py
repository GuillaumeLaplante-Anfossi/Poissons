'''
A quick script for generating diagonals using the shift method of SU from
"Comparing Diagonals on the Associahedra"


Use the two natural shift pairs
-  (L,R)    i.e the SU original
-  (L',R')  a symmetry of the prior pair, and we now know these are LA diagonals.
'''
from itertools import permutations, chain, combinations, product
from copy import deepcopy
import sys

def powerset(iterable):
    '''
    https://stackoverflow.com/questions/1482308/how-to-get-all-subsets-of-a-set-powerset
    returns a generator for the powerset
    '''
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def vertex_generator(n):
    '''
    returns a generator (iterable) for the permutations of size n,
    which correspond to vertices.
    '''
    return permutations(range(1,n+1))

def strong_complimentary_pair(v):
    '''
    returns the strong complimentary pair of a given vertex
    '''
    merge_desc = [[v[0]]]
    merge_asc = [[v[0]]]

    prev = v[0]
    for x in v[1:]:
        if x>prev:
            merge_asc[-1].append(x)
            merge_desc.append([x])
        if x<prev:
            merge_desc[-1].append(x)
            merge_asc.append([x])
        prev = x

    #Sort the subpartitions so in standard form
    merge_asc = [sorted(x) for x in merge_asc]
    merge_desc = [sorted(x) for x in merge_desc]

    return merge_desc, merge_asc

def R_shifts(lcp,j=0):
    '''
    returns all possibe R shifts of the left elem of complementary pair when j=0
    '''
    shifts = [lcp]

    if j == len(lcp)-1:
        #no right shift possible
        return shifts

    feasible_Ms = powerset(lcp[j][1:])

    for M in feasible_Ms:
        if len(M)==0:
            #If M is emptyset go to next.
            shifts = R_shifts(lcp,j+1)
        #min M_j > max A_{j+1}  then perform shift
        elif M[0]>lcp[j+1][-1]:
            shift = deepcopy(lcp) #copy the contents not the pointer
            shift[j] = list(set(lcp[j]) - set(M))
            shift[j+1]  = sorted(shift[j+1]+list(M))
            
            #Recurse
            shifts = shifts + R_shifts(shift,j+1)

    return shifts


def L_shifts(rcp,j=None):
    '''
    returns all possibe L shifts of the right elem of complementary pair

    We use the indexing B_1...B_q (instead of Bq...B1)
    Also refer to Ns as Ms
    '''
    shifts = [rcp]
    if j is None:
        j = len(rcp)-1
    if j==0:
        #No left shift possible
        return shifts

    feasible_Ms = powerset(rcp[j][1:])

    for M in feasible_Ms:
        if len(M)==0:
            #If M is emptyset go to next.
            shifts = L_shifts(rcp,j-1)
        #min M_j > max A_{j-1}  then perform shift
        elif M[0]>rcp[j-1][-1]:
            shift = deepcopy(rcp) #copy the contents not the pointer
            shift[j] = list(set(rcp[j]) - set(M))
            shift[j-1] = sorted(shift[j-1]+list(M))

            #Recurse
            shifts = shifts + L_shifts(shift,j-1)

    return shifts

def LR_shift_of_v(v):
    '''
    returns A_v x B_v
    '''
    lcp, rcp = strong_complimentary_pair(v)
    pairs = []
    for l in R_shifts(lcp):
        for r in L_shifts(rcp):
            pairs.append([l,r])

    return pairs

def diagonal_via_LR_shift(n):
    pairs = []
    for v in vertex_generator(n):
        pairs = pairs + LR_shift_of_v(v)
    return pairs


def R_prime_shifts(lcp,j=None):
    '''
    returns all possibe R shifts of the left elem of complementary pair when j=0
    '''
    shifts = [lcp]

    if j is None:
        j = len(lcp)-1

    if j == 0:
        #no right prime shift possible
        return shifts

    feasible_Ms = powerset(lcp[j][:-1])

    for M in feasible_Ms:
        if len(M)==0:
            #If M is emptyset go to next.
            shifts = R_prime_shifts(lcp,j-1)
        #max M_j < min A_{j-1} then perform shift
        elif M[-1]<lcp[j-1][0]:
            shift = deepcopy(lcp) #copy the contents not the pointer
            shift[j] = list(set(lcp[j]) - set(M))
            shift[j-1]  = sorted(shift[j-1]+list(M))
            
            #Recurse
            shifts = shifts + R_prime_shifts(shift,j-1)

    return shifts


def L_prime_shifts(rcp,j=None):
    '''
    returns all possibe L shifts of the right elem of complementary pair

    We use the indexing B_1...B_q (instead of Bq...B1)
    Also refer to Ns as Ms
    '''
    shifts = [rcp]
    if j is None:
        j = 0
    if j==len(rcp)-1:
        #No left shift possible
        return shifts

    feasible_Ms = powerset(rcp[j][:-1])

    for M in feasible_Ms:
        if len(M)==0:
            shifts = L_prime_shifts(rcp,j+1)
        #max M_j < min A_{j+1}  then perform shift
        elif M[-1]<rcp[j+1][0]:
            shift = deepcopy(rcp) #copy the contents not the pointer
            shift[j] = list(set(rcp[j]) - set(M))
            shift[j+1] = sorted(shift[j+1]+list(M))

            #Recurse
            shifts = shifts + L_prime_shifts(shift,j+1)

    return shifts

def LR_prime_shift_of_v(v):
    '''
    returns A_v x B_v
    '''
    lcp,rcp = strong_complimentary_pair(v)
    pairs = []
    for l in R_prime_shifts(lcp):
        for r in L_prime_shifts(rcp):
            pairs.append([l,r])

    return pairs

def diagonal_via_LR_prime_shift(n):
    pairs = []
    for v in vertex_generator(n):
        pairs = pairs + LR_prime_shift_of_v(v)
    return pairs

def SU_diag(n):
    return diagonal_via_LR_shift(n)

def LA_diag(n):
    return diagonal_via_LR_prime_shift(n)

if __name__ == "__main__":
    v = list(vertex_generator(3))[3]
    print(v)
    lcp,rcp = strong_complimentary_pair(v)
    print(lcp,rcp)
    print(R_prime_shifts([[2],[1,3,4]]))
    print()

    print(diagonal_via_LR_shift(3))
    print()
    print(diagonal_via_LR_prime_shift(3))

    for k in range(1,7):
        n1 = len(diagonal_via_LR_shift(k))
        n2 = len(diagonal_via_LR_prime_shift(k))
        print(n1,n2)


    #Printing out lists for comparison.
    for k in range(1,7):
        s = "./Examples/Shift/n={}.txt".format(k)
        with open(s, 'w') as f:
            sys.stdout = f # Change the standard output to file.
            print(diagonal_via_LR_shift(k))

        s = "./Examples/PrimeShift/n={}.txt".format(k)
        with open(s, 'w') as f:
            sys.stdout = f # Change the standard output to file.
            print(diagonal_via_LR_prime_shift(k))