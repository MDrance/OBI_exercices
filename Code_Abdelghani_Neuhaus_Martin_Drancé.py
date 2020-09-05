#Code de Abdelghani Neuhaus et Martin Drancé


################ IMPORTATION DES MODULES NECESSAIRES ################ 
import random

################ CREATION DES DONNEES ################ 
# Génération d'une séquence de longueur choisie
def generate_random_seq (num, alpha):                          
    seqtmp = []
    listalpha = list(alpha)
    for i in range(0, num):
        randnumber = random.randint(0,len(alpha)-1)
        seqtmp.append(listalpha[randnumber])
    seq = "".join(seqtmp)
    return seq

# Génération d'une liste de qualité (autant de qualités que de nucléotides dans le read)
def quality_values (read):
    qualitytmp = []
    for i in range (0, len(read)):
        proba = round(random.random(), 2)
        qualitytmp.append(str(proba))
    quality = " ".join(qualitytmp)
    return quality

# Ajout des adaptateurs aux extrémités des reads
def add_adaptateur (data_reads, adaptateurs):
    for i in data_reads:
        if data_reads.index(i) % 2 == 0:
            pos_to_keep = data_reads.index(i)
            data_reads.pop(data_reads.index(i))
            rdm = random.randint(0, 1)
            if rdm == 0:
                i = random.choice(adaptateurs) + i
                data_reads.insert(pos_to_keep, i)
            elif rdm == 1:
                rdm_reverse = random.choice(adaptateurs)[::-1]
                i = i + rdm_reverse
                data_reads.insert(pos_to_keep, i)
    return data_reads

# Extraction de string à des positions définies
def oneWord (seq, start, wlen):
    word = ""
    for i in range(len(seq)):
        if i >= start and i <= start + wlen-1:
            word += seq[i]
    return word 

# Génération des reads, associés à des qualités
def reads (seq):
    data_reads = []
    for i in range(len(seq)):
        if i == len(seq) - 20:
            break
        else:
            wlen = random.randint(20, 100)
            read = (oneWord (seq, i, wlen))
            qual = quality_values(read)
            data_reads.append(read)
            data_reads.append(qual)           
    return data_reads

# Suppresion des adaptateurs, qu'ils soient au début ou à la fin (reverse)
def FiltreAdaptateurs(readset, adaptateurs):
    for i in range(len(readset)):
        if i % 2 == 0:
            for a in range(len(adaptateurs)):
                if adaptateurs[a] in readset[i]:
                    value_to_keep = i
                    new_read = readset[i].replace(adaptateurs[a],"")
                    readset.pop(i)
                    readset.insert(value_to_keep, new_read)
                elif adaptateurs[a][::-1] in readset[i]:
                    value_to_keep = i
                    new_read = readset[i].replace(adaptateurs[a][::-1],"")
                    readset.pop(i)
                    readset.insert(value_to_keep, new_read)
    return readset

# Ecriture des données dans un fichier
def wfic(data_reads, file):
    ficw = open (file, "w")
    for data in data_reads:
        ficw.write (str(data))
        ficw.write ("\n")
    ficw.close()

# Lecture d'un fichier
def rfic (file):
    fic  = open(file, "r")
    line = fic.readline()
    data_reads = []
    while (line != ""):
        read = line.strip()
        data_reads.append(read)
        line = fic.readline()
        qual = line.split()
        qual_float = [float(i) for i in qual]
        data_reads.append(qual_float)
        line = fic.readline()
    fic.close()
    return data_reads

# Vérification des reads et de leur nature (ADN)
def isDNA (seq):
    for nuc in seq:
        if nuc != "a" and nuc != "t" and nuc != "c" and nuc != "g":
            if nuc != "A" and nuc != "T" and nuc != "C" and nuc != "G":
                return print ("Le read",seq, "n'est pas de l'ADN")

# Vérification de l'intervalle des qualités (entre 0 et 1)
def isquality (quality):
    for i in quality:
        if i < 0 or i > 1:
            return print ("Les valeurs de qualité", quality, "sont hors limites")

# Vérification de la longueur des reads (entre 20 et 100)
def isgoodlen (string):
    if len(string) < 20 or len(string) > 100:
        return print ("Le read", string, "n'a pas une longueur correcte")

# Vérification qu'un read et une liste de qualités associée aient la même longueur
def issamelen (read, qual):
    for i,j in zip(read, qual):
        if len(i) != len(j):
            return print (i, j, "n'ont pas la même taille")

# On check tout ce qui a été dit dans les 4 précédentes fonctions
def check_before_LoadReadSet(file):
    raw_data = rfic(file)
    read_liste = []
    qual_liste = []
    for i in range(len(raw_data)):
        if i%2 == 0:
            isDNA(raw_data[i])
            isgoodlen(raw_data[i])
            read_liste.append(raw_data[i])
        elif i%2 == 1:
            isquality(raw_data[i])
            qual_liste.append(raw_data[i])
    issamelen(read_liste, qual_liste)
    return read_liste, qual_liste


################ CHARGEMENT DES DONNEES PRECEDEMMENT TRAITEES ET SAUVEGARDEES ################ 
# Chargement des données après avoir vérifiés tous les éléments à prendre en compte (cf guide_fonctions.txt)
def LoadReadSet(file):
    reads = check_before_LoadReadSet(file)[0]
    qual = check_before_LoadReadSet(file)[1]
    read_set = {}
    read_id = 0
    for i,j in zip(reads, qual):
        read_set[read_id] = [i,j,len(i)]
        read_id += 1
    return read_set


################ FILTRATION DES DONNEES ################ 
# Filtrage des reads par les 2 extrémités selon la qualité de chaque nucléotide du read        
def FiltreExtremites(ReadSet, SeuilQualite):
    readset = {}
    read_list = []
    qual_list = []
    read_id = 0
    for key in ReadSet.keys():
        read_list.append(list(ReadSet[key][0]))
        qual_list.append(ReadSet[key][1])
    for read, qual in zip(read_list, qual_list):
        try:
            while qual[0] < SeuilQualite:
                read.pop(0)
                qual.pop(0) 
            while qual[-1] < SeuilQualite:
                read.pop(-1)
                qual.pop(-1) 
        except IndexError:
            read_list.pop(read_list.index(read))
            qual_list.pop(qual_list.index(qual))
    for i,j in zip(read_list, qual_list):
        readset[read_id] = ["".join(i),j,len(i)]
        read_id += 1
    return readset

# Filtrage des reads selon la qualité moyenne
def Filtre (ReadSet, Seuilqualite):
    readset = {}
    readlist = []
    quallist = []
    read_id = 0
    for key in ReadSet.keys():
        readlist.append(list(ReadSet[key][0]))
        quallist.append(ReadSet[key][1])
    for read, qual in zip(readlist, quallist):
        som = 0.0
        moy = 0.0
        for value in qual:
            som = som + value
            moy = som / len(qual)
        if moy < Seuilqualite:
            readlist.pop(readlist.index(read))
            quallist.pop(quallist.index(qual))
    for i,j in zip(readlist, quallist):
        readset[read_id] = ["".join(i),j,len(i)]
        read_id += 1
    return readset   

# Renvoie du reverse complément d'un read
def ReverseComplement(Read):
    table = Read.maketrans("ATCG", "TAGC")
    reverse = Read.translate(table)
    return reverse


################ RECONSTRUCTION DE LA SEQUENCE INITIALE ################ 
# Donne le chevauchement maximal entre 2 séquences
def lenchevauchementMax(seq1, seq2):
   m = [[0] * (1 + len(seq2)) for i in range(1 + len(seq1))]    
   count = 0                                              
   longueur_max = 0                                            
   for x in range(1, 1 + len(seq1)):                          
       for y in range(1, 1 + len(seq2)):                      
           if seq1[x - 1] == seq2[y - 1]:                      
               m[x][y] = m[x - 1][y - 1] + 1                
               if m[x][y] > count:
                   count = m[x][y]                          
                   longueur_max = x
           else:
               m[x][y] = 0                                  
   return (seq1[longueur_max-count:longueur_max])

# Donne la position du read qui a le chevauchement maximal dans un ReadSet par rapport à un read défini
def idchevauchementmax (ReadSet, seq):
    readlist = []
    for key in ReadSet.keys():
        readlist.append(ReadSet[key][0])
    best_id = 0
    best_len = 0
    for read in readlist:
        longueur = len(lenchevauchementMax(read, seq))
        key = readlist.index(read)
        if longueur > best_len:
            best_len = longueur
            best_id = key
    return print ("Meilleur chevauchement exacte avec le read",best_id,"avec",best_len,"chevauchements")

# Assemblage des séquences qui ont le meilleur chevauchement 
def fusionSequence(seq1, seq2): 
    new_seq = ""
    nuc_to_remove = lenchevauchementMax(seq1, seq2)
    seq2 = seq2.replace(nuc_to_remove,"")
    new_seq = seq1 + seq2
    return new_seq

# Assemblage de la séquence initale (ou du moins la plus proche possible)
def reconstruct(ReadSet, seq):
    final_seq = ""
    readlist = []
    for key in ReadSet.keys():
        readlist.append(ReadSet[key][0])
    if not final_seq:
        final_seq = fusionSequence(readlist[0], readlist[1])
    if final_seq:
        for i in range(2, len(readlist)-1):
            final_seq = fusionSequence(final_seq,readlist[i])
            if len(final_seq) > len(seq):
                break
    return final_seq

# Affichage du pourcentage de similarité entre la séquence initiale et la séquence reconstruite après traitement
def accuracy(seq_init, final_seq):
    count = 0
    for i in range(len(seq_init) - 1):
        if seq_init[i] == final_seq[i]:
            count += 1
    return count/len(seq_init)



################################################
################### MAIN #######################

adaptateurs = ["TATATATA", "AAAAAT", "CCCCC"]                       #On genère les adaptateurs que l'on utilisera                                                             
seq = generate_random_seq(500, "ATCG")                              #On genere la séquence de base aléatoirement avec une longueur définie et les nucléotides ATCG                         
reads = reads(seq)                                                  #On creer les reads à partir de notre sequence et on leur attribue des qualités avec quality_values()    
# print ("Les reads sont", reads,"\n")                                #                   
adapt = add_adaptateur(reads, adaptateurs)                          #On ajoute un adaptateur à chacun de nos reads, au début ou à la fin                               
# print ("reads avec adapt", adapt,"\n")                              #                            
nonadapt = FiltreAdaptateurs(adapt, adaptateurs )                   #On enlève les adaptateurs de nos reads, du début ou de la fin                                       
# print ("reads sans adapt", nonadapt,"\n")                           #                                
wfic(nonadapt, "nonadapt.txt")                                      #On enregistre notre ensemble read/qualité dans un fichier text                    
readset = LoadReadSet("nonadapt.txt")                               #On lis nos reads et leur qualités depuis le fichier .txt, on test qu'ils sont valides avec check_before_loadset()                            
# print ("read set après vérifications", readset,"\n")                #                                            
dicofinal = FiltreExtremites(readset, 0.3)                          #On met dans un dictionnaire nos reads, leurs qualités et un ID dans un dictionnaire après filtrage des extremités                                
# print ("readset avec extremité filtre", dicofinal,"\n")             #                                            
readfinal = Filtre(dicofinal, 0.2)                                  #On filtre les reads dont la qualité moyenne est trop basse                       
# print ("readset avec reads mauvaise qualite filtre", readfinal)     #                                                    
idchev = idchevauchementmax(readfinal, "ATCTCGATCGTTAGCT")          #On trouve parmis nos reads celui qui a le meilleur chevauvement avec celui passé en argument                                           
recons = reconstruct(readfinal, seq)                                #On reconstruit notre séquence de base à partir des reads filtrés                           
print ("séquence initiale", seq,"\n")                               #        
print ("séquence reconstruite", recons)                             #    
accu = accuracy(seq, recons)                                        #On compare notre ADN initial et celui reconstruit                   
print ("La précision est de", accu)                                                          





