from time import process_time
import sys

#adaptation au cas permutations seulement et ménage de printemps
n=4
# j'ai supprimé min, max qui existent à priori déjà

ll=Permutations(n) # liste des permutations sous forme de partition ordonnée

# paires I,J de même cardinal avec min I < min J
def ijsubsets(n,k):
    l = []
    for x in Subsets(n,2*k):
        x2 = x.symmetric_difference(Set([min(x)]))
        for y in Subsets(x2,k):
            x3 = x.symmetric_difference(y)
            l.append([sorted(x3),sorted(y)])
    return l

vv=[a for k in range(1,n//2+1) for a in ijsubsets(n,k)] 
# liste des couples (I,J) d'ensembles de taille k 

## fonction auxiliaire (procédure) qui effectue une permutation circulaire sur les éléments d'indices pairs d'une liste an "avançant les éléments vers la gauche"

def rotEven(l):
    aux=l[0];
    lastInd=0
    for i in range(2,len(l),2):
        l[i-2]=l[i]
        lastInd=i
    l[lastInd]=aux
    
## fonction qui étant données deux permutations p1 et p2 de la forme partitions ordonnées et une paire de I,J vérifiant les conditions de Guillaume détermine si la condition poisson est vérifiée.
## renvoie un booléen

def poissonB(p1,p2,IJ):
    p1red=[] #intersection de p1 avec I \cup J
    p2red=[] #intersection de p2 avec I \cup J
    nextSet=1 #permet de s'assurer de l'alternance voulue pour la condition poisson. Vaut 0 quand le prochain élément à voir est dans I et 1 quand il est dans J
    IUJ=IJ[0]+IJ[1]
    for i in range(len(p1)):
        if p1[i] in IJ[nextSet]:
            p1red+=[p1[i]]
            nextSet= (nextSet+1)%2 # on a trouvé un élément de l'un : on change pour trouver un élément de l'autre
        elif p1[i] in IJ[(nextSet+1)%2]:
            return False # si on rencontre un élément qui n'est pas du bon ensemble, on s'arrête aussitôt de calculer
        if p2[i] in IUJ:
            p2red+=[p2[i]] # on continue de calculer l'intersection de p2 avec IUJ
    if len(p1red)!=len(p2red) or len(p1red)!=len(IUJ):
        return False
    else:
        i1=p1red[1] # on prend le min
        if not(i1==min(IUJ)):
            return False
        p1Circ=[p1red[i] for i in range(-1,len(p1red)-1)] # on effectue la transformation circulaire de p1red
        rotEven(p1Circ) # on décale les indices des éléments de I
        rotEven(p1Circ) # et encore une fois
        return p2red==p1Circ


# fonction qui teste si le couple de permutations p1/p2 vérifie la contre-condition de Guillaume (c'est à dire n'est pas dans la diagonale pour I,J en témoin). Les permutations sont supposées de même longueur
def testGuillaume(p1,p2,I,J):
    nbJmIp1=0 # différence du nombre d'éléments de J et du nombre d'éléments de I vus dans p1 à l'instant t
    nbImJp2=0 # différence du nombre d'éléments de I et du nombre d'éléments de J vus dans p2 à l'instant t
    for i in range(len(p1)):
        if p1[i] in J:
            nbJmIp1 += 1 
        elif p1[i] in I:
            nbJmIp1 += -1
        if p2[i] in I:
            nbImJp2 += 1 
        elif p2[i] in J:
            nbImJp2 += -1
        if nbJmIp1<0 or nbImJp2<0:
            return False
    return True
    
# fonction qui teste si deux partttions ordonénes sont dans la diagonale
def testDiag(p1,p2,I,J):
    nbImJp1=0 # différence du nombre d'éléments de J et du nombre d'éléments de I vus dans p1 à l'instant t
    nbJmIp2=0 # différence du nombre d'éléments de I et du nombre d'éléments de J vus dans p2 à l'instant t
    for i in range(len(p1)):
        for j in p1[i]:
            if j in I:
                nbImJp1 += 1 
            elif j in J:
                nbImJp1 += -1
        if nbImJp1>0:
            return True
    for i in range(len(p2)):
        for k in p2[i]:
            if k in J:
                nbJmIp2 += 1 
            elif k in I:
                nbJmIp2 += -1
        if nbJmIp2>0:
            return True
    return False
    
# test direct de partitions ordonnées
def test(p1,p2):
    for [I,J] in vv: 
        if not(testDiag(p1,p2,I,J)):
            return False
    return True
    
## fonction qui vérifie le critère poisson sur toutes les permutations
## affiche toutes les paires de permutations comparables pour l'ordre faible, mais pas dans la diagonale
#vv : liste des (I,J)
#ll : liste des permutations
#n : n courant
##n=5, affiche 90
## le test avec n=6 et ll=permutations met 140sec, affiche 3312
## le test avec n=7 et ll=ordre faible met quelques heures, affiche 122 400
## après benchmark, affiche
## pour n= 7
## longueur de vv 175
## ll et vv initialisés
## longueur de ll 672697
## ... décompte jusque 647630 --> à quoi correspond-il ?
## la condition poisson est bien vérifiée partout !
## 16.345780801000018 250.77963907799997 267.125419879
## pour n=8, le calcul s'interrompt au moment de l'initialisation de la liste ll
## /usr/share/sagemath/bin/sage-python : ligne 2 : 46739 Processus arrêté      sage -python "$@"

def poissMain(n): 
    print("pour n=",n)
    start=process_time()
    good=1
    vv=[a for k in range(2,n//2+1) for a in ijsubsets(n,k)] 
    print("longueur de vv",len(vv))
    ll=[i for i in Permutations(n).weak_lattice().relations()]
    print("ll et vv initialisés")
    print("longueur de ll",len(ll))
    middle=process_time()
    count=0
    for [x,y] in ll:
        if x!=y:
            np=0 # indique le nombre d'éléments sur lesquels un test poisson qui échoue
            IJ=[] # paire des IJ qui font échouer
            for [I,J] in vv: #on teste sur tous les (I,J)
                if testGuillaume(x,y,I,J) :
                    good=0 # good vaut 1 quand x et y sont dans la diagonale
                    IJ+=[[I,J]]; # ajout de la paire problématique
            if good==0 : # échoue le test, mais pas pour les inversions
                poissonV=[]
                poissonF=[]
                for ij in IJ:
                    if poissonB(x,y,[ij[0], ij[1]]):
                        poissonV+=[[Set(ij[0]), Set(ij[1])]]
                    else:
                        poissonF+=[[Set(ij[0]), Set(ij[1])]]
                if len(poissonF)!=0:
                    np=len(poissonF)
                    for ij in poissonF:
                        abloc=False
                        for ab in poissonV:
                            if ab[0].issubset(ij[0]) and ab[1].issubset(ij[1]):
                                abloc=True
                        if abloc:
                            np-=1
                        else:
                            print(x,y, poissonV, poissonF, ij)
                    if np!=0:
                        print("il y a un hic")
                count=count+1
                if count%1000==0:
                    print(count)
    print (count)
    if np==0:
        print("la condition poisson est bien vérifiée partout !")
    end=process_time()
    print(middle-start,end-middle, end-start)



# fonction qui génère tous les couples qu'il nous faut        
def ordreFaible(n):
    l_perm=[]
    p=Permutations(n)
    for x in p:
        for y in x.permutohedron_greater():
#            if y!=x:
# décommenter quand on cherche les motifs poissons
            l_perm+=[[x,y]]
    return l_perm        
        
                
## fonction comme la précédente, mais avec une autre manière de générer ll
#vv : liste des (I,J)
#ll : liste des permutations
#n : n courant
# si la conjecture est vraie, on doit retrouver le même nombre que précédemment
# un peu moins efficace que la version Main sur la création de paires
## pour n= 4
## longueur de vv 3
## longueur de ll 127
## nbre de couples qui ne sont pas dans la diagonale : 2
## à comparer avec 2
## la condition poisson est bien vérifiée partout !
## 0.001861265999999695 0.000726762000000214 0.002588027999999909
## pour n=5
## longueur de vv 15
## longueur de ll 1779
## nbre de couples qui ne sont pas dans la diagonale : 90
## à comparer avec 90
## la condition poisson est bien vérifiée partout !
## 0.031438687999999715 0.05062205200000047 0.08206074000000019
## pour n= 6
##longueur de vv 55
##longueur de ll 30991
## nbre de couples qui ne sont pas dans la diagonale : 3312
## à comparer avec 3312
## la condition poisson est bien vérifiée partout !
## 1.0201610000003711 3.5473550360002264 4.5675160360005975
## pour n=7
## longueur de vv 175
## longueur de ll 667657
## nbre de couples qui ne sont pas dans la diagonale : 122400
## à comparer avec 122400
## la condition poisson est bien vérifiée partout !
## 73.30185978400004 263.410195403 336.71205518700003
## pour n=8
## longueur de vv 525
## longueur de ll 17 511 003
## ...
## 17511003
## la condition poisson est bien vérifiée partout !
## 8506.674139542 21494.009500223 30000.683639764997

def poissFaible(n): 
    with open('poissonFFseul'+str(n)+'.txt', 'w') as f:
        sys.stdout=f
        print("pour n=",n)
        start=process_time()
        vv=[a for k in range(2,n//2+1) for a in ijsubsets(n,k)]
        print("longueur de vv",len(vv))
        ll=ordreFaible(n)
        print("longueur de ll",len(ll))
        middle=process_time()
        good=1
        count=0
        for [x,y] in ll:
            good=1
            np=0 # indique le nombre d'éléments sur lesquels un test poisson qui échoue
            IJ=[] # paire des IJ qui font échouer
            for [I,J] in vv: #on teste sur tous les (I,J)
                if testGuillaume(x,y,I,J) :
                    good=0 # good vaut 1 quand x et y sont dans la diagonale
                    IJ+=[[I,J]]; # ajout de la paire problématique
#        print("x=", x,", y= ",y,", IJ=", IJ)
            if good==0 : # échoue le test, mais pas pour les inversions
                poissonV=[]
                poissonF=[]
                for ij in IJ:
                    if poissonB(x,y,[ij[0], ij[1]]):
                        poissonV+=[[Set(ij[0]), Set(ij[1])]]
                    else:
                        poissonF+=[[Set(ij[0]), Set(ij[1])]]
                if len(poissonF)!=0:
                    np=len(poissonF)
                    for ij in poissonF:
        	            abloc=False
        	            for ab in poissonV:
        	                if ab[0].issubset(ij[0]) and ab[1].issubset(ij[1]):
        	                    abloc=True
        	            if abloc:
        	                np-=1
        	            else:
        	                print(x,y, poissonV, poissonF, ij)
                    if np!=0:
                        print("il y a un hic")
                    print("x=",x," , y=",y,", poissonV=",poissonV,", poissonF=", poissonF)
#            print("="*10)
                count=count+1
##                print("x=", x, ", y=", y, ", poissonV=", poissonV, "poissonF=", poissonF)
        print ("nbre de couples qui ne sont pas dans la diagonale :",count)
        diag=[1, 1, 3, 17, 149, 1809, 28399, 550297, 12732873, 343231361, 10576764251, 367054970721, 14173669352413, 602974492511377, 28027436035348359, 1413479599558432169, 76879014760731439889, 4486205132570631391617, 279595430611791210216883] # OEIS A213507 : nombre de couples de permutations dans la diagonale
        wo=[1, 1, 3, 17, 151, 1899, 31711, 672697, 17551323, 549500451, 20246665349, 864261579999, 42190730051687, 2329965898878307] # OEIS A007767 : nombre d'intervalles de l'ordre faible
        print("à comparer avec",wo[n]-diag[n])
        if np==0:
            print("la condition poisson est bien vérifiée partout !")
        end=process_time()
        print(middle-start,end-middle, end-start)

# fonction qui renvoie la liste des [I,J] vérifiant la conjecture d'Hugues
# fonction cassée à réparer... (test n=4)

def IJreduit(n):
    vvDeb=[a for k in range(2,n//2+1) for a in ijsubsets(n,k)] 
    nbx= 0 # nombre de x vu dans l'intervalle [1,t] où t est la position courante
    nby= 0 # nombre de y vu dans l'intervalle [1,t] où t est la position courante
    vvFin=[]
    for [x,y] in vvDeb:
        nbx= 0 # nombre de x vu dans l'intervalle [1,t] où t est la position courante
        nby= 0 # nombre de y vu dans l'intervalle [1,t] où t est la position courante
        cond=True
        for t in range(min(x)+1, min(y)):
            if t in x:
                cond=False
                break
        if cond:
            for t in range(min(y),n+1):
                if t in x :
                    nbx+=1
                if t in y:
                    nby+=1
                if nbx>=nby:
                    cond=False
                    break
        if cond:
            vvFin+=[[x,y]] # si la condition d'Hugues est  vérifié, on l'ajoute à la liste
    return vvFin
            

# procédure qui vérifie si la conjecture est vérifiée avec un argument de dénombrement
# sage:poissOracle(3)                          
# n= 3
# On veut trouver  0  paires
# longueur de vv 0
# longueur de ll 11
# nbre de couples qui vérifie poisson : 0
# 0.0008505190000960283 0.00012186299977656745 0.0009723819998725958
# sage:poissOracle(4)             
# n= 4
# On veut trouver  2  paires
# longueur de vv 1
# longueur de ll 127
# nbre de couples qui vérifie poisson : 2
## nombre de paires (I,J) utilisées :  1
## [[[1, 4], [2, 3]]]
## 0.0018400100000235398 0.00047271399989767815 0.002312723999921218
# sage: poissOracle(5)           
# n= 5
# On veut trouver  90  paires
## longueur de vv 5
## longueur de ll 1779
## nbre de couples qui vérifie poisson : 90
## nombre de paires (I,J) utilisées :  5
## [[[2, 5], [3, 4]], [[1, 5], [2, 4]], [[1, 4], [2, 3]], [[1, 5], [2, 3]], [[1, 5], [3, 4]]]
## 0.034213825000051656 0.025045007000016994 0.05925883200006865
# sage: poissOracle(6)                                  
# n= 6
# On veut trouver  3312  paires
## longueur de vv 17
## longueur de ll 30991
## nbre de couples qui vérifie poisson : 3312
## nombre de paires (I,J) utilisées :  17
## [[[3, 6], [4, 5]], [[2, 6], [3, 5]], [[2, 5], [3, 4]], [[2, 6], [3, 4]], [[2, 6], [4, 5]], [[1, 6], [2, 5]], [[1, 5], [2, 4]], [[1, 6], [2, 4]], [[1, 4], [2, 3]], [[1, 4, 6], [2, 3, 5]], [[1, 5], [2, 3]], [[1, 6], [2, 3]], [[1, 5, 6], [2, 3, 4]], [[1, 6], [3, 5]], [[1, 5], [3, 4]], [[1, 6], [3, 4]], [[1, 6], [4, 5]]]
## 0.9965422919999583 1.4687693739999759 2.465311665999934
# sage: poissOracle(7)                                  
# n= 7
# On veut trouver  122400  paires
# longueur de vv 49
# longueur de ll 667657
# trouvés : 100000
# nbre de couples qui vérifie poisson : 122400
# nombre de paires (I,J) utilisées :  49
# [[[4, 7], [5, 6]], [[3, 7], [4, 6]], [[3, 6], [4, 5]], [[3, 7], [4, 5]], [[3, 7], [5, 6]], [[2, 7], [3, 6]], [[2, 6], [3, 5]], [[2, 7], [3, 5]], [[2, 5], [3, 4]], [[2, 5, 7], [3, 4, 6]], [[2, 6], [3, 4]], [[2, 7], [3, 4]], [[2, 6, 7], [3, 4, 5]], [[2, 7], [4, 6]], [[2, 6], [4, 5]], [[2, 7], [4, 5]], [[2, 7], [5, 6]], [[1, 7], [2, 6]], [[1, 6], [2, 5]], [[1, 7], [2, 5]], [[1, 5], [2, 4]], [[1, 5, 7], [2, 4, 6]], [[1, 6], [2, 4]], [[1, 7], [2, 4]], [[1, 4], [2, 3]], [[1, 4, 7], [2, 3, 6]], [[1, 4, 6], [2, 3, 5]], [[1, 4, 7], [2, 3, 5]], [[1, 5], [2, 3]], [[1, 5, 7], [2, 3, 6]], [[1, 6], [2, 3]], [[1, 7], [2, 3]], [[1, 6, 7], [2, 4, 5]], [[1, 6, 7], [2, 3, 5]], [[1, 5, 6], [2, 3, 4]], [[1, 5, 7], [2, 3, 4]], [[1, 6, 7], [2, 3, 4]], [[1, 7], [3, 6]], [[1, 6], [3, 5]], [[1, 7], [3, 5]], [[1, 5], [3, 4]], [[1, 5, 7], [3, 4, 6]], [[1, 6], [3, 4]], [[1, 7], [3, 4]], [[1, 6, 7], [3, 4, 5]], [[1, 7], [4, 6]], [[1, 6], [4, 5]], [[1, 7], [4, 5]], [[1, 7], [5, 6]]]
# 79.87540812600014 119.49560694699994 199.37101507300008
# n= 8 (lancé à 21h43)
# On veut trouver  4818450  paires
# longueur de vv 131
# longueur de ll 17511003
# trouvés : 100000
# ....trouvés : 4800000
# nbre de couples qui vérifie poisson : 4818450
# nombre de paires (I,J) utilisées :  131
# [[[5, 8], [6, 7]], [[4, 8], [5, 7]], [[4, 7], [5, 6]], [[4, 8], [5, 6]], [[4, 8], [6, 7]], [[3, 8], [4, 7]], [[3, 7], [4, 6]], [[3, 8], [4, 6]], [[3, 6], [4, 5]], [[3, 6, 8], [4, 5, 7]], [[3, 7], [4, 5]], [[3, 8], [4, 5]], [[3, 7, 8], [4, 5, 6]], [[3, 8], [5, 7]], [[3, 7], [5, 6]], [[3, 8], [5, 6]], [[3, 8], [6, 7]], [[2, 8], [3, 7]], [[2, 7], [3, 6]], [[2, 8], [3, 6]], [[2, 6], [3, 5]], [[2, 6, 8], [3, 5, 7]], [[2, 7], [3, 5]], [[2, 8], [3, 5]], [[2, 5], [3, 4]], [[2, 5, 8], [3, 4, 7]], [[2, 5, 7], [3, 4, 6]], [[2, 5, 8], [3, 4, 6]], [[2, 6], [3, 4]], [[2, 6, 8], [3, 4, 7]], [[2, 7], [3, 4]], [[2, 8], [3, 4]], [[2, 7, 8], [3, 5, 6]], [[2, 7, 8], [3, 4, 6]], [[2, 6, 7], [3, 4, 5]], [[2, 6, 8], [3, 4, 5]], [[2, 7, 8], [3, 4, 5]], [[2, 8], [4, 7]], [[2, 7], [4, 6]], [[2, 8], [4, 6]], [[2, 6], [4, 5]], [[2, 6, 8], [4, 5, 7]], [[2, 7], [4, 5]], [[2, 8], [4, 5]], [[2, 7, 8], [4, 5, 6]], [[2, 8], [5, 7]], [[2, 7], [5, 6]], [[2, 8], [5, 6]], [[2, 8], [6, 7]], [[1, 8], [2, 7]], [[1, 7], [2, 6]], [[1, 8], [2, 6]], [[1, 6], [2, 5]], [[1, 6, 8], [2, 5, 7]], [[1, 7], [2, 5]], [[1, 8], [2, 5]], [[1, 5], [2, 4]], [[1, 5, 8], [2, 4, 7]], [[1, 5, 7], [2, 4, 6]], [[1, 5, 8], [2, 4, 6]], [[1, 6], [2, 4]], [[1, 6, 8], [2, 4, 7]], [[1, 7], [2, 4]], [[1, 8], [2, 4]], [[1, 4], [2, 3]], [[1, 4, 8], [2, 3, 7]], [[1, 4, 7], [2, 3, 6]], [[1, 4, 8], [2, 3, 6]], [[1, 4, 6], [2, 3, 5]], [[1, 4, 6, 8], [2, 3, 5, 7]], [[1, 4, 7], [2, 3, 5]], [[1, 4, 8], [2, 3, 5]], [[1, 5], [2, 3]], [[1, 5, 8], [2, 3, 7]], [[1, 5, 7], [2, 3, 6]], [[1, 5, 8], [2, 3, 6]], [[1, 6], [2, 3]], [[1, 6, 8], [2, 3, 7]], [[1, 7], [2, 3]], [[1, 8], [2, 3]], [[1, 7, 8], [2, 5, 6]], [[1, 7, 8], [2, 4, 6]], [[1, 4, 7, 8], [2, 3, 5, 6]], [[1, 7, 8], [2, 3, 6]], [[1, 6, 7], [2, 4, 5]], [[1, 6, 8], [2, 4, 5]], [[1, 6, 7], [2, 3, 5]], [[1, 6, 8], [2, 3, 5]], [[1, 7, 8], [2, 4, 5]], [[1, 7, 8], [2, 3, 5]], [[1, 5, 6], [2, 3, 4]], [[1, 5, 6, 8], [2, 3, 4, 7]], [[1, 5, 7], [2, 3, 4]], [[1, 5, 8], [2, 3, 4]], [[1, 5, 7, 8], [2, 3, 4, 6]], [[1, 6, 7], [2, 3, 4]], [[1, 6, 8], [2, 3, 4]], [[1, 7, 8], [2, 3, 4]], [[1, 6, 7, 8], [2, 3, 4, 5]], [[1, 8], [3, 7]], [[1, 7], [3, 6]], [[1, 8], [3, 6]], [[1, 6], [3, 5]], [[1, 6, 8], [3, 5, 7]], [[1, 7], [3, 5]], [[1, 8], [3, 5]], [[1, 5], [3, 4]], [[1, 5, 8], [3, 4, 7]], [[1, 5, 7], [3, 4, 6]], [[1, 5, 8], [3, 4, 6]], [[1, 6], [3, 4]], [[1, 6, 8], [3, 4, 7]], [[1, 7], [3, 4]], [[1, 8], [3, 4]], [[1, 7, 8], [3, 5, 6]], [[1, 7, 8], [3, 4, 6]], [[1, 6, 7], [3, 4, 5]], [[1, 6, 8], [3, 4, 5]], [[1, 7, 8], [3, 4, 5]], [[1, 8], [4, 7]], [[1, 7], [4, 6]], [[1, 8], [4, 6]], [[1, 6], [4, 5]], [[1, 6, 8], [4, 5, 7]], [[1, 7], [4, 5]], [[1, 8], [4, 5]], [[1, 7, 8], [4, 5, 6]], [[1, 8], [5, 7]], [[1, 7], [5, 6]], [[1, 8], [5, 6]], [[1, 8], [6, 7]]]
# 8680.734895293 6312.357728191999 14993.092623484998
##(précedemment :  8506.674139542 21494.009500223 30000.683639764997)
# n= 10
# On veut trouver  9669901098  paires
# longueur de vv 869
# arrêté le calcul de ll au bout de 2h...


def poissOracle(n):
    print("n=",n)
    diag=[1, 1, 3, 17, 149, 1809, 28399, 550297, 12732873, 343231361, 10576764251, 367054970721, 14173669352413, 602974492511377, 28027436035348359, 1413479599558432169, 76879014760731439889, 4486205132570631391617, 279595430611791210216883] # OEIS A213507 : nombre de couples de permutations dans la diagonale
    wo=[1, 1, 3, 17, 151, 1899, 31711, 672697, 17551323, 549500451, 20246665349, 864261579999, 42190730051687, 2329965898878307] # OEIS A007767 : nombre d'intervalles de l'ordre faible
    print("On veut trouver ",wo[n]-diag[n]," paires")
    start=process_time()
    vv=IJreduit(n) 
    print("longueur de vv",len(vv))
    ll=ordreFaible(n)
    print("longueur de ll",len(ll))
    middle=process_time()
    good=1
    count=0
    IJ=[]
    for [x,y] in ll:
        good=1
        for [I,J] in vv: #on teste sur tous les (I,J)
            if poissonB(x,y,[I,J]):
                good=0 # good vaut 0 quand x et y vérifie la condition de poisson forte
                if [I,J] not in IJ:
                    IJ+=[[I,J]]
                break
        if good==0:
            count=count+1
            if count%100000==0:
                print("trouvés :",count)
    print ("nbre de couples qui vérifie poisson :",count)
    end=process_time()
    print("nombre de paires (I,J) utilisées : ",len(IJ))
    print(IJ)
    print(middle-start,end-middle, end-start)

## fonction qui renvoie les sous-ensembles de [I,J] vérifiant la condition de Hugues



def goodSubset(I,J):
    aux=[[A,B] for k in range(2,len(I)) for A in Combinations(I,k).list() for B in Combinations(J,k).list() ]
    nba=0
    nbb=0
    res=[]
    for [A,B] in aux:
        if min(A)>min(B):
            cond=False
        else:
            nba= 0 # nombre de x vu dans l'intervalle [1,t] où t est la position courante
            nbb= 0 # nombre de y vu dans l'intervalle [1,t] où t est la position courante
            cond=True
            for t in range(min(A)+1, min(B)):
                if t in A:
                    cond=False
                    break
            if cond:
                for t in range(min(B),n+1):
                    if t in A :
                        nba+=1
                    if t in B:
                        nbb+=1
                    if nba>=nbb:
                        cond=False
                        break
        if cond:
            res+=[[A,B]] # si la condition d'Hugues est  vérifié, on l'ajoute à la liste
    return res
    

## Dernier essai d'accélération : utiliser le fait que si (I,J) bug pour la condition de Guillaume, alors un de ses sous-ensembles doit vérifier la condition poisson
## bug --> à voir !

def poissTGV(n): 
    print("pour n=",n)
    diag=[1, 1, 3, 17, 149, 1809, 28399, 550297, 12732873, 343231361, 10576764251, 367054970721, 14173669352413, 602974492511377, 28027436035348359, 1413479599558432169, 76879014760731439889, 4486205132570631391617, 279595430611791210216883] # OEIS A213507 : nombre de couples de permutations dans la diagonale
    wo=[1, 1, 3, 17, 151, 1899, 31711, 672697, 17551323, 549500451, 20246665349, 864261579999, 42190730051687, 2329965898878307] # OEIS A007767 : nombre d'intervalles de l'ordre faible
    print("On veut trouver ",wo[n]-diag[n]," paires")
    start=process_time()
    vv=[a for a in ijsubsets(n,n//2)] 
    print("longueur de vv",len(vv))
    ll=ordreFaible(n)
    print("longueur de ll",len(ll))
    middle=process_time()
    good=True  # vaut True si la paire est dans la diagonale
    count=0 
    poisson=False # a-t-on trouvé un poisson
    for [x,y] in ll:
        good=True
        poisson=False  
        for [I,J] in vv: #on teste sur tous les (I,J)
            if testGuillaume(x,y,I,J) : 
                count+=1
#                print("x=",x,", y=",y,", [I,J]=",[I,J])
                good=False # good vaut 1 quand x et y sont dans la diagonale
                if poissonB(x,y,[I,J]):
                    poisson=True
##                    print("poisson avec ",[I,J])
                else:
##                    print("pas de poisson avec [I,J]")
                    for AB in goodSubset(I,J):
                        if poissonB(x,y,AB):
                            poisson=True
##                            print("poisson avec ",AB)
                            break
                if not(poisson):
                    print("cas sans sous-poisson trouvé pour x=",x," et y=",y, "[I,J]=",I,J)
                break
    print ("nbre de couples qui ne sont pas dans la diagonale :",count)
    end=process_time()
    print(middle-start,end-middle, end-start)

## fonction qui compte le nombre de couples dans la diagonale en fonction de leurs dimensions

def nbDiag(n):
    print("n=",n)
    diag=[1, 1, 3, 17, 149, 1809, 28399, 550297, 12732873, 343231361, 10576764251, 367054970721, 14173669352413, 602974492511377, 28027436035348359, 1413479599558432169, 76879014760731439889, 4486205132570631391617, 279595430611791210216883] # OEIS A213507 : nombre de couples de permutations dans la diagonale
    print("On veut trouver ",diag[n]," paires")
    start=process_time()
    vv=IJreduit(n) 
#    print("longueur de vv",len(vv))
    ll=ordreFaible(n)
#    print(ll)
    middle=process_time()
    good=1
    count=0
    for [x,y] in ll:
#        print(x,y)
        good=1
        for [I,J] in vv: #on teste sur tous les (I,J)
            if poissonB(x,y,[I,J]):
                good=0 # good vaut 0 quand x et y vérifie la condition de poisson forte
                break
        if good==1:
            count=count+1
            if count%100000==0:
                print("trouvés :",count)
#        print(good)
    print ("nbre de couples danns la diagonale :",count)
    end=process_time()
    print(middle-start,end-middle, end-start)

