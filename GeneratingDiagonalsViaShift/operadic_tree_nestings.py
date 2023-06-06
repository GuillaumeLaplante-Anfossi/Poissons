import sys
from itertools import permutations, chain, combinations, product
from copy import deepcopy

from diagonals_via_shift import SU_diag, LA_diag

class Tree:
    def __init__(self, data = None, children = []):
        self.data = data
        self.children = children

    #https://stackoverflow.com/questions/20242479/printing-a-tree-data-structure-in-python
    def __str__(self, level=0):
        ret = "\t"*level+repr(self.data)+"\n"
        for child in self.children:
            ret += child.__str__(level+1)
        return ret

    def __repr__(self):
        return '<tree node representation>'

def non_empty_powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1,len(s)+1))

def increment_tree(t,v):
    '''Update the data of all vertices of t to data += v.'''
    base = t
    t.data = t.data + v
    for child in t.children:
        increment_tree(child,v)

    return base

def PT(n):
    '''The set of planar trees with n non-root vertices labelled under infix order.

    Exploiting catalan recurrence for generation.
        C = 1 + z C^2 
    '''
    if n == 0:
        return [Tree(0)]

    trees = []
    for k in range(n):
        rec_1 = PT(k)
        rec_2 = PT(n-1-k)
        for t1 in rec_1:
            for t2 in rec_2:
                t = deepcopy(t1)
                t.children.append(increment_tree(deepcopy(t2),k+1))
                trees.append(t)

    return trees

def connected_to_root(t):
    '''Returns a set of all subtrees of t which are connected to root.'''
    if t.children == []:
        return [] #If no children then have no sub-trees
    #Recurse
    connected_to_children = []
    for child in t.children:
        child_and_edge_to_root = [[child.data]] #the edge to child is a trivial subtree
        for x in connected_to_root(child):
            x.append(child.data)
            child_and_edge_to_root.append(x) #the edge to child, and each subtree of child is a subtree
        connected_to_children.append(child_and_edge_to_root)

    #Combine
    all_conn_to_root = []
    #Choose a non-emptyset of children to be connected to
    for choose_conn_children in non_empty_powerset(connected_to_children):
        #Construct the product of all sub-trees from this choice.
        for p in product(*choose_conn_children):
            sub = []
            for e in p:
                sub += e
            all_conn_to_root.append(sub)  

    return all_conn_to_root

def compute_nests(t):
    '''Computes all nests of a tree.
    This is equivalent to computing all connected subtrees of t, 
    which is equivalent to the union of all vertices in t,
        subtrees of t connected to root 
    '''
    #base case
    nests = connected_to_root(t)
    #recursion
    for child in t.children:
        sub_nests = compute_nests(child)
        if sub_nests != []:
            nests += sub_nests
    return nests

def check_nesting_of_tree(d,t):
    '''
    Checks that sigma, tau is a nesting of tree t.

    By construction meets compatibility, and trivial nest.
    So just need to check each element is a nest.
    '''
    sigma,tau = d
    N = compute_nests(t)
    N = [sorted(x) for x in N]
    N.sort()

    #Manipulate into sigma_1 | sigma_1 \cup sigma_2|...
    sigma_nests = []
    for i in range(1,len(sigma)+1):
        sigma_nests.append(sorted(list(chain(*sigma[:i]))))

    #check
    for block in sigma_nests:
        if block not in N:
            return False

    #repeat
    tau_nests = []
    for i in range(1,len(tau)+1):
        tau_nests.append(sorted(list(chain(*tau[:i]))))
    for block in tau_nests:
        if block not in N:
            return False

    return True

def proj_nesting(sigma,t):
    '''
    Given an ordered partition sigma, projects it onto a tree t.

    Ex tree from email:

    '''
    proj = []
    N = compute_nests(t)
    N = [sorted(x) for x in N]
    N.sort()

    print(N)

    #Manipulate into sigma_1 | sigma_1 \cup sigma_2|...
    merged = []
    for i in range(1,len(sigma)+1):
        merged.append(sorted(list(chain(*sigma[:i]))))

    for block in merged:
        if block in N:
            proj.append(block)

def construct_tree_nesting_dict(n, diag = "SU"):
    if diag == "SU":
        computed_diag = SU_diag(n)
    elif diag == "LA":
        computed_diag = LA_diag(n)

    pt = PT(n)
    tree_nesting_dict = {}
    for t in pt:
        #Use strings to uniquely identify trees
        s = str(t)
        tree_nesting_dict[s] = [] 
    #Incredibly weird bug, can only iterate through computed diag once...
    #Maybe did something weired with iterators, can't remember.
    #Fine as is, but finicky.
    for d in computed_diag:
        sigma,tau = d
        for t in pt:      
            if (not check_proj_nesting(sigma,t)) or (not check_proj_nesting(tau,t)):
                continue
            s = str(t)
            tree_nesting_dict[s].append(d)

    return tree_nesting_dict

def check_proj_nesting(sigma,t):
    ''''''
    return len(sigma) == len(proj_nesting(sigma,t))

def proj_nesting(sigma,t):
    '''
    Given an ordered partition sigma, projects it onto a tree t.
    Ex tree from email:
    '''
    proj = []
    N = compute_nests(t)
    N = [sorted(x) for x in N]
    N.sort()

    #Manipulate into sigma_1 | sigma_1 \cup sigma_2|...
    merged = []
    for i in range(1,len(sigma)+1):
        merged.append(sorted(list(chain(*sigma[:i]))))

    #print(N)
    #print(sigma)
    #print(merged)
    #print("loop")

    #Break each merged block into minimal nests
    min_nesting = []
    for block in merged:
        #print("block",block)
        min_nests_block = []

        start = 0
        for i in range(1,len(block)+1):
            #print(block[start:i],block[start:i] in N)
            if block[start:i] in N:
                continue
            else:
                #print("found", i, block[start:(i-1)])
                min_nests_block.append(block[start:(i-1)])
                start=i-1
        #Wrap up
        if start<len(block)+1:
            min_nests_block.append(block[start:])
        #inefficient
        for m in min_nests_block:
            if m not in min_nesting:
                min_nesting.append(m)
    #print("min_nesting",min_nesting)
    #print("\n")
    return min_nesting

if __name__ == '__main__':

    orig_stdout = sys.stdout
    f = open('bijection_check.txt', 'w')
    sys.stdout = f


    n = 5
    pt = PT(n)
    '''
    t = pt[0]
    print(t)
    print(proj_nesting([[1, 4], [2, 3]],t))
    print(check_proj_nesting([[1, 4], [2, 3]],t))
    print(proj_nesting([[1, 4], [2], [3]],t))
    print(check_proj_nesting([[1, 4], [2], [3]],t))
    print(proj_nesting([[3, 4], [2], [1]],t))
    print(check_proj_nesting([[3, 4], [2], [1]],t))

    print(proj_nesting([[1, 3], [2], [4]],t))
    print(check_proj_nesting([[1, 3], [2], [4]],t))
    
    exit()
    

    nesting_dict_SU = construct_tree_nesting_dict(n, diag = "SU")
    nesting_dict_LA = construct_tree_nesting_dict(n, diag = "LA")
    t = pt[33]
    s = str(t)
    print(s)
    print("countLA",len(nesting_dict_LA[s]))
    print("countSU",len(nesting_dict_SU[s]))
    diffLA = []
    diffSU = []
    for x in nesting_dict_LA[s]:
        if x not in nesting_dict_SU[s]:
            diffLA.append(x)
    for x in nesting_dict_SU[s]:
        if x not in nesting_dict_LA[s]:
            diffSU.append(x)
    print("in LA but not in SU",len(diffLA))
    print("in SU but not in LA",len(diffSU))

    sortedLA = [[sorted(x[0]),sorted(x[1])] for x in nesting_dict_LA[s]]
    sortedSU = [[sorted(x[0]),sorted(x[1])] for x in nesting_dict_SU[s]]

    print("In SU and reordering of blocks is not in LA")
    for i in range(len(nesting_dict_SU[s])):
        if sortedSU[i] not in sortedLA:
            print(nesting_dict_SU[s][i])
    print("In LA and reordering of blocks is not in SU")
    for i in range(len(nesting_dict_LA[s])):
        if sortedLA[i] not in sortedSU:
            print(nesting_dict_LA[s][i])

    print("\n")
    print("All LA")
    for x in nesting_dict_LA[s]:
        print(x)

    print("\n")
    print("All SU")
    for x in nesting_dict_LA[s]:
        print(x)

    exit()
    '''

    for n in range(3,5+1):
        print("For n=",n)
        nesting_dict_SU = construct_tree_nesting_dict(n, diag = "SU")
        nesting_dict_LA = construct_tree_nesting_dict(n, diag = "LA")
        trees = PT(n)
        for i in range(len(trees)):
            t = trees[i]
            s = str(t)

            if len(nesting_dict_SU[s]) != len(nesting_dict_LA[s]):
                print("Tree:")
                print(s)
                #print(compute_nests(t))
                print("SU Count:",len(nesting_dict_SU[s]))
                print("LA Count:",len(nesting_dict_LA[s]))

                diffLA = []
                diffSU = []
                for x in nesting_dict_LA[s]:
                    if x not in nesting_dict_SU[s]:
                        diffLA.append(x)
                for x in nesting_dict_SU[s]:
                    if x not in nesting_dict_LA[s]:
                        diffSU.append(x)

                print("disjoint",len(diffLA),len(diffSU))

            else:
                print("For tree i=",i+1,"shared count=",len(nesting_dict_SU[s]))
            

    sys.stdout = orig_stdout
    f.close()