README.TXT de Abdelghani Neuhaus et Martin Drancé


On importe le module random qui va nous servir à générer des valeurs aléatoires comprises dans des intervalles définis et aussi à choisir aléatoirement entre A, T, C et G lors de la création des reads.

D'abord, nous présenterons les fonctions une par une, puis nous expliquerons où elles ont été utilisées et appelées dans le programme principal.




################## Fonction generate_random_seq(num,alpha) ##################

Cette fonction sert à générer une séquence aléatoire de longueur 'num' et composée de lettres mis dans un string 'alpha'. Par exemple, si num=100 et alpha='ABCD', une séquence de 100 lettres tirées aléatoirement à chaque fois entre A, B, C et D sera générée. La fonction renvoie la séquence créée.




################## Fonction quality_values(read) ##################

Cette fonction va générer une suite de nombres à virgule de longueurs 'read', c'est-à-dire de la longueur du string placé en argument. Ces nombres seront tirés au sort aléatoirement dans un intervalle [0 , 1]. Par exemple, si read = 'hello', la fonction renvoie [0.87 0.11 0.72 0.22 0.99]




################## Fonction add_adaptateur(data_reads,adaptateurs) ##################

Cette fonction renvoie une liste des reads générés précédement et y ajoute un adaptateur soit au début, soit à la fin. Si l'adaptateur est ajouté à la fin, il est à l'envers.




################## Fonction oneWord(seq, start, wlen) ##################

Cette fonction extrait un mot 'seq' de longueur 'wlen' à la position 'start'.




################## Fonction reads(seq) ##################

Cette fonction prend en argument un string 'seq'. Elle renvoie une liste contenant des sous-string (ce sont les reads) issus de 'seq', de longueur comprise entre 20 et 100 (chaque read est généré avec la fonction oneWord) et une liste de qualité associée, générée avec la fonction quality_values.




################## Fonction FiltreAdaptateurs(readset,adaptateurs) ##################

Cette fonction prend en argument un readset (donc reads et les qualités associées) et la liste des adaptateurs. La fonction permet de parcourir chaque read et de regarder la présence de chaque adaptateurs au début ou à la fin de celui-ci. Lorsque l'adaptateur est trouvé sur un read, il est supprimé. Il a été aussi pris en compte le cas où l'on a le brin reverse (donc l'inverse du motif de l'adaptateur).




################## Fonction wfic(data_reads,file) ##################

Cette fonction écrit sur un fichier file.txt les données que l'on a. Les éléments sont écrits ligne par ligne (une ligne pour un read, puis une ligne pour les qualités, puis une ligne pour le read suivant, etc...)




################## Fonction rfic(file) ##################

Cette fonction lit et récupère les données d'un fichier nommé file.txt.




################## Fonction isDNA(seq) ##################

Cette fonction prend en argument une séquence 'seq' et vérifie que les lettres la composant soient bien A, T, C ou G (minuscule ou majuscule). Si ce n'est pas le cas, un message d'avertissement est renvoyé. Sinon, rien n'est envoyé.




################## Fonction isquality(quality) ##################

Cette fonction vérifie que les qualités contenues dans les différentes listes soient bien des nombres compris dans l'intervalle [0,1]. Si ce n'est pas le cas, un message d'avertissement est affiché.




################## Fonction isgoodlen(string) ##################

Cette fonction vérifie que la longueur du read est comprise dans l'intervalle dans lequel ils ont été générés (ici, il faut qu'ils aient une longueur comprise entre 20 et 100).




################## Fonction issamelen(read,qual) ##################

Cette fonction vérifie que la longueur du read et la longueur de la liste de qualité associée soient les mêmes.




################## Fonction check_before_LoadReadSet(file) ##################

Cette fonction fait appel aux 4 précédentes fonctions et va lire un fichier de type TXT et checker successivement si les reads sont bien de l'ADN, si leur longueur est comprise entre 20 et 100 nucléotides, si les qualités associées sont des valeurs comprises entre 0 et 1 et si les reads et les liste de qualité associées ont la même longueur. Si tout est bon, on renvoit une liste des reads et une liste des qualités.




################## Fonction LoadReadSet(file) ##################

Cette fonction fait appel à la fonction précédente. Elle va simplement mettre dans un dictionnaire un read et la liste de qualité associée. L'identifiant du read sera donc la clé du dictionnaire.




################## Fonction FiltreExtremites(ReadSet, SeuilQualite) ##################

Cette fonction prend en argument un jeu de données sous forme de dictionnaire et un seuil de qualité. Pour chaque read, on parcourt chaque nucléotide et la qualité associée. Si cette dernière est inférieure au seuil entré en paramètre, la qualité et le nucléotide associée sont supprimés. Le filtrage s'arrête une fois qu'une des qualité est supérieure au seuil. Il est aussi possible que tout le read soit filtré. Le filtrage commence aux extrémités du read. Il est ainsi renvoyé un dictionnaire mis-à-jour.




################## Fonction Filtre(ReadSet, SeuilQualite) ##################

Cette fonction prend en argument un jeu de données sous forme de dictionnaire et un seuil de qualité. Pour chaque read, on parcourt la liste de qualité associée. On fait la moyenne des qualités composant la liste. Si cette moyenne est inférieure au seuil, le read et la liste associée sont supprimés. On renvoit un dictionnaire mis-à-jour.




################## Fonction ReverseComplement(Read) ##################

Cette fonction prend en argument une séquence ADN et renvoie son reverse complement




################## Fonction chevauchementMax(ReadSet,Read) ##################

Cette fonction prend en argument un read et le ReadSet. On va parcourir chaque read du ReadSet et renvoie l’identifiant du read du ReadSet ayant le plus grand chevauchement avec le Read passé en paramètre de la fonction.




################## Fonction ReverseComplement(Read) ##################

Cette fonction prend en argument une séquence ADN et réalise l’assemblage de reads obtenus après les étapes de filtrage. Elle renvoie une séquence unique correspondant à la séquence initiale (plus ou moins)




################## Fonction lenchevauchementMax(seq1,seq2) ##################

Cette fonction prend en argument 2 séquences et revoie la longueur du chevauchement maximum. Cette fonction crée un tableau à double entrée avec en horizontal la sequence 1 et à la verticale la séquence 2. Chaque case est identifiable par des coordonnées (comme une matrice). 




################## Fonction idchevauchementmax(ReadSet,seq) ##################

Cette fonction prend en argument le ReadSet et le read d'intérêt. On parcourt le ReadSet et on fait appel à chaque itération à la fonction précédente pour renvoyer la longueur du chevauchement max entre le read d'intérêt et les reads du ReadSet. On retient la position et la meilleure longueur. C'est ce qui est renvoyé par la fonction.




################## Fonction fusionSequence(seq1,seq2) ##################

Cette fonction prend en argument 2 reads. On appelle la fonction lenchevauchemenMax pour trouver le chevauchement maximal entre les 2 séquences. Ensuite, on retire ce motif (le chevauchement) de la 2è séquence et on ajoute ce qui reste à la première séquence.
Par exemple : on a 'personne' et 'nenuphar' : le plus long chevauchement est 'ne'. On le retire de nenuphar, ce qui donne 'nuphar'. On ajoute ceci à 'personne' : on obtient 'personnenuphar'. C'est ce qui est renvoyé par la fonction.




################## Fonction reconstruct(ReadSeq,seq) ##################

Cette fonction prend en argument le ReadSet et un read. Elle fait appel à la fonction précédente et réalise la fusion de séquence à chaque tour jusqu'à que tous les reads aient été assemblés. On renvoit la séquence finale, qui est plus ou moins proche de la séquence initiale générée.




################## Fonction accuracy(seq_init,final_seq) ##################

Cette fonction prend en argument la séquence initiale et la séquence générée par la fonction précédente. On regarde le pourcentage d'identité nucléotide par nucléotide.


Nous avons remarqué que le pourcentage est souvent supérieure à 95%, mais des fois, il est de 25 % environ. Cela correspond au fait qu'il y'a un seul décalage dans la séquence. Ainsi, les séquences sont identiques à un nucléotide près. Mais, étant donné qu'ici on regarde une similarité nucléotide par nucléotide, il n'y a pas homologie.

Par exemple : 
'iloveparis'
'loveparis'

Ici, on aurait une homologie de 0%, alors qu'il y'a qu'un nucléotide de décalage. Il faut donc prendre en compte ceci lors de la lecture du pourcentage d'homologie renvoyée par la fonction.








########################################################################################################
########################################################################################################
########################################################################################################

Programme principal :


1) d'abord, on crée des adaptateurs (longueurs entre 6 et 10 nucléotides)

2) ensuite, on génère une séquence composée de A, T, C et G de 500 nucléotides

3) on fait appel à la fonction qui génère les reads à partir de la séquence précédemment générée et on leur associe des qualités

4) on ajoute les adaptateurs au début ou à la fin, dans le sens normal ou reverse.

5) on supprime les adaptateurs, peu importe leurs sens

6) on sauvegarde nos data dans un fichier TXT

7) on récupère ces data et on check les différents éléments présentés dans les fonctions précédemment (reads composés de A, T, C ou G ; longueur de reads correspond à celles des qualités ; les qualités sont comprises entre 0 et 1 ; les reads ont une longueur comprise entre 20 et 100 nucléotides)

8) on met les reads dans un dictionnaire avec leur qualité. L'identifiant correspond à la clé du dictionnaire.

9) on filtre les reads par les extrémités jusqu'à que les nucléotides aux extrémités aient une qualité supérieure au seuil

10) on filtre les reads une deuxième fois, cette fois par rapport à leur qualité moyenne

11) on regarde quel read du ReadSet a le meilleur chevauchement avec le read placé en argument

12) on reconstruit une séquence entière à partir des différents reads

13) on regarde le pourcentage d'identité avec la séquence initale


















