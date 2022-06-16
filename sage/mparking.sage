
def desc(l):
    d = []
    for i in range(0,len(l)-1):
        if l[i]>l[i+1]:
            d.append(i)
    return d


# variante des fonctions de parking: suite de nb strictement positifs telle 
# qu'une fois triée on reste  sous 234...(n+1).
def park(n):
    l = []
    for i in range( n , sum(range(n+1,1,-1))+1 ):
        for p in Partitions(i,length=n,outer=range(n+1,1,-1)):
            print(p)
            l += Permutations(p)
    return l


# add( q**len(desc(list(p)+[1])) for p in park(2) )
# c'était un essai mais ça ne donne pas les nombres qu'on cherche

# le multigraphe Kn / e
def gr(n):
    G = graphs.CompleteGraph(n)
    G.allow_multiple_edges(true)
    G.contract_edge(G.edges()[0])
    return G

# gr(4).spanning_trees_count()

