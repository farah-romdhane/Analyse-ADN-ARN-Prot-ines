#TP1
# Date : 25-03-2024
# Nom : 
# Farah ROMDHANE (Matricule :20288662) 
# Sara EL AHOUEL(Matricule : 20257092)


#Fonction antisens: Cette fonction prend en entrée un brin d'ADN et renvoie son brin complémentaire antisens.

def antisens(brinAdn):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'} # Dictionnaire pour la correspondance des nucléotides  
    complementaireadn = ""  
    for i in brinAdn:  
        complementaireadn += complement[i]  
    return complementaireadn[::-1]  # Retourner le brin complémentaire inversé pour le lire de gauche à droite
    
#Fonction trouveDebut : Recherche tous les codons de départ sur un brin d’ADN et renvoie un tableau contenant les positions du premier nucléotide de chacun des codons

def trouveDebut(brinAdn):
    debutPositions = []
    for i in range (len(brinAdn)-2) :
        if brinAdn [i:i+3]== "TAC" : # Vérifier si la séquence courante est un codon de début
            debutPositions.append(i)   
    return debutPositions #renvoie la liste des positions de début


#Fonction trouveFin: Même chose que la fonction précédente mais renvoie un tableau avec les positions de tous les codons de terminaison

def trouveFin(brinAdn):
    finPositions = []
    # Parcours du brin d'ADN
    for i in range(len(brinAdn)-2):
        # Vérification si un des codons de fin est présent à cette position
        if brinAdn[i:i+3] in ['ATT', 'ATC', 'ACT']:
            # Ajout de la position du premier nucléotide du codon de fin
            finPositions.append(i)
    return finPositions

"""Fonction trouveGene : Prend en paramètre un tableau contenant les positions de tous les codons de départ et un autre tableau contenant 
les positions de tous les codons de terminaison pour un brin d’ADN et renvoie un tableau de tuples contenant la liste des gènes 
(début et fin) trouvés sur un brin."""
    
def trouveGene(debut, fin):
    genes = []                      
    for debutGene in debut:          
        for finGene in fin:            # Pour chaque codon de début, chercher un codon de fin correspondant
            if finGene > debutGene and (finGene - debutGene) % 3 == 0 :
                genes.append((debutGene, finGene))     # Si oui, ajouter le tuple (debutGene, finGene) à la liste des gènes identifiés
                # Après avoir trouvé le premier codon de fin valide pour ce codon de début, passer au prochain codon de début
                break
    return genes

#Fonction extraireSequenceADN :Extraie les séquences de gènes à partir d'un brin d'ADN basé sur les positions de début et de fin des gènes 
    
def extraireSequenceADN(brinAdn, gene):
    sequencesGenes = []
    for debutGene, finGene in gene:
        segmentAdnDuGene = brinAdn[debutGene:finGene]   # Extraire chaque segment de gène à partir du brin d'ADN
        sequencesGenes.append(segmentAdnDuGene)
    return sequencesGenes

#Fonction transcrire : Transcrit des séquences de gènes d'ADN en ARN.

def transcrire(sequencesGenes):
    correspondanceArn = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
    sequencesArn = []  # Liste pour stocker les séquences d'ARN transcrit
    for sequenceGene in sequencesGenes:# Parcourir chaque séquence de gène dans la liste
        sequenceArn = ""  
        for nucl in sequenceGene:  # Parcourir la séquence de gène actuelle
            if nucl in correspondanceArn:
                sequenceArn += correspondanceArn[nucl]  # Transcription
        sequencesArn.append(sequenceArn)  
    return sequencesArn

#Fonction traduire : Prend en paramètre un brin d’ARN (chaine de caractères) et affiche la protéine sous forme d’une chaine de caractères

def traduire(sequenceArn):
    codons_aa = {
        "UUU": "Phénylalanine",
        "UUC": "Phénylalanine",
        "UUA": "Leucine",
        "UUG": "Leucine", 
        "CUU": "Leucine", 
        "CUC": "Leucine",
        "CUA": "Leucine",
        "CUG": "Leucine",
        "AUU": "Isoleucine",
        "AUC": "Isoleucine",
        "AUA": "Isoleucine", 
        "AUG": "Méthionine ",
        "GUU": "Valine",
        "GUC": "Valine",
        "GUA": "Valine", 
        "GUG": "Valine", 
        "UCU": "Sérine",
        "UCC": "Sérine",
        "UCA": "Sérine",
        "UCG": "Sérine",
        "CCU": "Proline",
        "CCC": "Proline",
        "CCA": "Proline", 
        "CCG": "Proline",
        "ACU": "Thrénine",
        "ACC": "Thrénine",
        "ACA": "Thrénine",
        "ACG": "Thrénine",
        "GCU": "Alanine", 
        "GCC": "Alanine",
        "GCA": "Alanine",
        "GCG": "Alanine", 
        "UAU": "Tyrosine",
        "UAC": "Tyrosine",
        "UAA": "Stop", 
        "UAG": "Stop",
        "CAU": "Histidine",
        "CAC": "Histidine",
        "CAA": "Glutamine",
        "CAG": "Glutamine",
        "AAU": "Asparagine",
        "AAC": "Asparagine",
        "AAA": "Lysine", 
        "AAG": "Lysine",
        "GAU": "Aspartate",
        "GAC": "Aspartate",
        "GAA": "Glutamate",
        "GAG": "Glutamate",
        "UGU": "Cystéine",
        "UGC": "Cystéine",
        "UGA": "Stop",
        "UGG": "Tryptophane",
        "CGU": "Arginine",
        "CGC": "Arginine",
        "CGA": "Arginine", 
        "CGG": "Arginine",
        "AGU": "Sérine",
        "AGC": "Sérine",
        "AGA": "Arginine",
        "AGG": "Arginine",
        "GGU": "Glycine",
        "GGC": "Glycine",
        "GGA": "Glycine",
        "GGG": "Glycine"
}
    proteines = []  # Liste pour stocker les séquences de protéines
    for sequenceArn in sequencesArnTranscrit:
        proteine = ["Méthionine (Start)"]  # Liste pour stocker les acides aminés de la protéine en cours de traduction
        for i in range(3, len(sequenceArn), 3):
            codon = sequenceArn[i:i+3]   # Extraire le codon courant.
            if codon in codons_aa:
                acideAmine = codons_aa[codon]
                if acideAmine == "Stop":
                    break             # Arrêter la traduction si le codon est un codon stop.
                proteine.append(acideAmine)
        proteines.append("-".join(proteine))# Ajouter la chaîne  traduite à la liste des protéines, séparés par des -
    return proteines

#Fonction affiche_proteine:traduit des séquences d'ARN en séquences de protéines et les affiche sous forme des lettres séparées par des tirets.

def afficheProteine(sequencesArnTranscrit):
    lettreAa = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}
    proteinesLettre = []       # Liste pour stocker toutes les séquences de protéines
    for sequenceArn in sequencesArnTranscrit:
        proteineLettre = []    # Liste pour stocker les aa de la protéine en cours de traduction
        for i in range(0, len(sequenceArn), 3):
            codon = sequenceArn[i:i+3]
            
            if codon in lettreAa:
                acideAmine = lettreAa[codon]
                if acideAmine == "*":       # Arrêter la traduction si un codon stop est trouvé
                    break
                proteineLettre.append(acideAmine)
        proteinesLettre.append("-".join(proteineLettre))  
    return proteinesLettre  # Retourner la liste des protéines traduites

#Affichage des protéines à l’aide de la tortue:

import turtle
clear(340,600)

#Fonction carre: Dessine un carré sur une grille basée sur un indice donné.

def carre(longueurCote, indice):
    colonne = indice % 15
    ligne = indice // 15
    #positionner le carré et le dessine à cette position
    x = colonne * longueurCote - 150  # Ajustez selon votre préférence
    y = -ligne * longueurCote + 270   # Ajustez selon votre préférence
    pu()
    goto(x, y)
    pd()
    for _ in range(4):    #dessiné le carré
        fd(longueurCote)
        rt(90)
    # se positioner au centre du carré pour ecrire la lettre 
    pu()
    goto(x + longueurCote / 2, y - longueurCote/ 2 )
    pd()


"""afficher_proteines: Prend une liste de séquences d'acides aminés et la longueur du côté de chaque carré comme paramètres. Pour chaque
acide aminé dans chaque séquence, la fonction trouve la position appropriée pour dessiner le carré, utilise la fonction carre pour le dessiner,
et écrit la lettre de l'acide aminé au centre. """

def afficherProteines(séquenceLetaa, longueurCote):
    nombre = 0                                  # L'indice du carré à dessiner
    for sequence in séquenceLetaa:
        for lettre in sequence.replace('-', ''):  #enlever le - entre les lettres
            carre(longueurCote, nombre)
            write(lettre)                          # écrire la lettre
            nombre += 1
            if nombre % 15 == 0 and nombre != 0 :  # Nouvelle ligne après 15 carrés
                nombre += -nombre % 15
        nombre += -nombre % 15 + 15 if nombre % 15 != 0 else 15  # Sauter une ligne entière à la fin de chaque protéine
        
#Programme Principal 
#appel des fonctions

# Brin d'ADN fourni

brinAdn="TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGACG\
GCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGAGTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGAATGCCAGCCAGC\
CAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGA\
ACTCGACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCATCCCAGCGATACCC\
AGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAG\
CCAGCGAACTCGTCTGCGTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGC\
GATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGTATGCCAGCC\
AGCATCCCAGCGA"

brinComplementaire=antisens(brinAdn)

debut1=trouveDebut(brinAdn) #les debuts qui se trouvent sur le brinAdn

fin1=trouveFin(brinAdn) #les fins qui se trouvent sur le brinAdn

geneadn =trouveGene(debut1, fin1) #les genes qui se trouvent sur le brinAdn

debut2=trouveDebut(brinComplementaire)#les debuts qui se trouvent sur le brin Complementaire

fin2=trouveFin(brinComplementaire)#les fins qui se trouvent sur le brin Complementaire

genecomplementaire=trouveGene(debut2, fin2)#les genes qui se trouvent sur le brin Complementaire

sequencesGenesAdn =extraireSequenceADN(brinAdn, geneadn) #sequence d'adn des genes pour le BrinAdn 

sequencesGenesComplementaire =extraireSequenceADN(brinComplementaire, genecomplementaire)#sequence d'adn des genes pour le Brin Complementaire

sequencesArnTranscrit = transcrire(sequencesGenesComplementaire)#la transcription 

proteine = traduire(sequencesArnTranscrit)#la traduction 

str = ' '.join(proteine)#convertir la liste en une chaine de caractéres

print("les proteines sur le brin complementaire sont : \n" ,str)

proteinLetaa = afficheProteine(sequencesArnTranscrit)

p1 = proteinLetaa 

sequencesArnTranscrit = transcrire(sequencesGenesAdn)

proteine = traduire(sequencesArnTranscrit)

str = ' '.join(proteine)

print("les proteines sur le brin adn  sont : \n",str)

proteinLetaa = afficheProteine(transcrire(sequencesGenesAdn))

p2=proteinLetaa

p=p2+p1 #concatenation des deux listes 

afficherProteines(p, 20)


#Tests unitaires
def testAntisens ():
    assert antisens ("ACGTTG") == "CAACGT"
    assert antisens ("TACTACGGTC") == "GACCGTAGTA"
    assert antisens ("GCGCGCGC") == "GCGCGCGC"
testAntisens()  

def testTrouveDebut():
    assert trouveDebut("AATTTACTTGGG") == [4]
    assert trouveDebut("AATGGCTTACTTGGG") == [7]
    assert trouveDebut("TACTTGGG") == [0]
testTrouveDebut()

def testTrouveFin():
    assert trouveFin("AATTTACTTGGGATC") == [1,5,12]
    assert trouveFin("AATGGCTTACTTGGGATC") == [8,15]
    assert trouveFin("TACTTGGG") == [1]
testTrouveFin()

def testTrouveGene():
    assert trouveGene([2,67],[9,15,100]) == [(67,100)]
    assert trouveGene([45,113,406],[15,135,293]) == [(45,135),(113,293)]
    assert trouveGene([25,690],[6,15,100,130]) == [(25,100)]
testTrouveGene()

def  testTranscrire():
    assert transcrire("TACTTCCCCCGGATT") == ['A', 'U', 'G', 'A', 'A', 'G', 'G', 'G', 'G', 'G', 'C', 'C', 'U', 'A', 'A']
    assert transcrire ("TACGGTTGCCTGACT")== ['A', 'U', 'G', 'C', 'C', 'A', 'A', 'C', 'G', 'G', 'A', 'C', 'U', 'G', 'A']
    assert transcrire ("TACGGTTACACT") == ['A', 'U', 'G', 'C', 'C', 'A', 'A', 'U', 'G', 'U', 'G', 'A']
testTranscrire()

def testTraduire():
    assert traduire(["GACCGU"]) == ["Méthionine (Start)-Arginine"]
    assert traduire(["AUGAUGGGU"]) == ["Méthionine (Start)-Méthionine-Glycine"]
    assert traduire(["AUGCUUGUAGGU"]) == ["Méthionine (Start)-Leucine-Leucine-Valine-Glycine"]
#testTraduire()