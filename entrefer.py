#!/usr/bin/env python3
import csv
import numpy as np
from scipy import integrate as it
from scipy.integrate import simps
from scipy import signal
from math import *
import matplotlib.pyplot as plt



def lecture_texte(nom_fichier):
    """
    en cours
    """


def lecture_csv():
    """
    lit le fichier .csv de haut en bas et de gauche à droite et renvoie le nombre
    de lignes, le nombre de colonnes, le numero de la colonne contenant le temps,
    le numero de la colonne contenant le toptour, les numeros des colonnes
    contenant des fems, les numeros des colonnes contenant des entrefers
    """
    id_temps = 0
    id_fem = 0
    id_TT = 0
    id_entrefer = 0
    les_fems = []
    les_entrefers = []
    nb_lignes = 0
    file = open("Auzat_new.csv","r")
    test = csv.reader(file)
    for row in test:
        nb_colonnes = (len(row))
        nb_lignes += 1
    file = open("Auzat_new.csv","r")
    ligne = file.readline()
    for j in range(nb_colonnes): # tant que element n'est pas egal a la derniere colonne
        element = ligne.split(',')[j]
        colonne_temps = element.find("(s)")
        colonne_fem = element.find("(V)")
        colonne_TT = element.find("TT")
        colonne_entrefer = element.find("(mm)")
        if colonne_temps != -1:
            id_temps = j
        if colonne_fem != -1:
            id_fem = j
            les_fems.append(id_fem)
        if colonne_TT != -1:
            id_TT = j
        if colonne_entrefer != -1:
            id_entrefer = j
            les_entrefers.append(id_entrefer)
    file.close()
    return(nb_colonnes, nb_lignes, id_temps, les_fems, id_TT, les_entrefers)



def analyse_temps(nb_lignes, id_temps):
    """
    lit la colonne du temps, replace le point de depart à 0s et convertit en ms
    pour plus de lisibilité et renvoie une liste avec tous les instants t
    """
    t = []
    file=open("Auzat_new.csv","r")
    test=csv.reader(file)
    ligne = file.readline()
    # on lit la 1ere ligne
    # on lit les autres lignes de données
    for i in range(2, nb_lignes + 1):
        ligne = file.readline()
        # on s'occupe d'abord du temps
        element = ligne.split(',')[id_temps]
        if i == 2:
            origine_des_temps = float(element)
        num = float(element)
        remis = (num - origine_des_temps)*1000
        # on se place par rapport à
        # l'origine du temps et on convertit en ms
        t.append(remis)
    file.close()
    return t



def analyse_TT(nb_lignes, id_TT):
    toptour = []
    file=open("Auzat_new.csv","r")
    test=csv.reader(file)
    ligne = file.readline()
    # on lit la 1ere ligne
    # on lit les autres lignes de données
    for i in range(2, nb_lignes + 1):
        ligne = file.readline()
        # pour le top tour
        element3 = ligne.split(',')[id_TT]
        tt = float(element3)
        toptour.append(tt)
    file.close()
    return toptour



def analyse_entrefer_k(nb_lignes, les_entrefers, k):
    """
    lit la k-ième colonne qui contient des infos sur un capteur d'entrefer
    et replace toutes ces valeurs dans une liste entrefer_k
    """
    entrefer_k = []
    file=open("Auzat_new.csv","r")
    test=csv.reader(file)
    ligne = file.readline()
    # on lit la 1ere ligne
    # on lit les autres lignes de données
    for i in range(2, nb_lignes + 1):
        ligne = file.readline()
        # pour la fem k
        element2 = ligne.split(',')[k]
        fe2 = float(element2)
        entrefer_k.append(fe2)
    file.close()
    return entrefer_k



def calcul_seuil(toptour, les_t):
    """
    calcule le seuil à partir duquel on déclenche le toptour et renvoie
    ce seuil et le sens du toptour (vers le haut ou vers le bas)
    """
    somme = 0
    dessus = []
    dessous = []
    ecart_max = 0
    sens = 0
    for i in toptour:
        somme += i
    moyenne = somme / (len(toptour))
    for j in toptour:
        ecart = abs(j - moyenne)
        if ecart > ecart_max:
            ecart_max = ecart
            id_max = j
        if j > moyenne:
            dessus.append(j)
        elif j < moyenne:
            dessous.append(j)
    if len(dessus) > len(dessous):
        seuil = id_max + (ecart_max * 3 / 3)  # on peut remplacer le 1er 3 par un 2
        sens = 1
    elif len(dessus) < len(dessous):
        seuil = id_max - (ecart_max * 3 / 3)
        sens = -1
    return (seuil, sens)



def tour_machine(toptour, les_t, seuil, sens):
    """
    calcule les instants pendant lesquels on franchit le seuil du toptour,
    retourne une liste avec les instants t qui correspondent à un
    tour machine
    """
    clef = 0
    dernier_tt = 0
    les_tt = []
    if sens == 1:
        for i in range(len(toptour)):
            if toptour[i] <= seuil and clef == 0:
                les_tt.append(i)
                dernier_tt = les_t[i]
                clef = 1
            ecart = les_t[i] - dernier_tt
            if ecart >= 10:
                clef = 0
    if sens == -1:
        for i in range(len(toptour)):
            if toptour[i] >= seuil and clef ==0:
                les_tt.append(i)
                dernier_tt = les_t[i]
                clef = 1
            ecart = les_t[i] - dernier_tt
            if ecart >= 10:
                clef = 0
    return(les_tt)



def selection_entrefer_k(entrefer_k, les_t, les_tt, seuil, sens, N):
    nb_pics = N * (len(les_tt) + 2)
    les_entrefers_k = [[] for _ in range(nb_pics + 1)]
    premier_tic = les_tt[0]
    dernier_tic = les_tt[-1]
    id_entrefer = -1
    reset = 0
    if sens == 1:
        """
        le cas qu'on a le plus souvent
        """
        for i in range(premier_tic, dernier_tic + 1):
            if entrefer_k[i] >= seuil and reset == 0:
                reset = 1
                id_entrefer += 1
            if entrefer_k[i] >= seuil and reset == 2:
                reset = 0
            if entrefer_k[i] <= seuil and reset != 0:
                reset = 2
                les_entrefers_k[id_entrefer].append(entrefer_k[i])
    if sens == -1:
        """
        l'autre sens
        """
        for i in range(premier_tic, dernier_tic + 1):
            if entrefer_k[i] <= seuil and reset == 0:
                reset = 1
                id_entrefer += 1
            if entrefer_k[i] <= seuil and reset == 2:
                reset = 0
            if entrefer_k[i] >= seuil and reset != 0:
                reset = 2
                les_entrefers_k[id_entrefer].append(entrefer_k[i])
    les_sommets_k = []
    for t in les_entrefers_k:
        if len(t) != 0:
            les_sommets_k.append(t)
    return(les_sommets_k)



def calcul_valeur_moyenne(liste, ecart_entre_t):
    ecart_gauche = 1.5624
    ecart_droit = 1.01556
    nb_gauche = round(ecart_gauche / ecart_entre_t)
    nb_droit = round(ecart_droit / ecart_entre_t)
    new_list = liste[nb_gauche:-nb_droit]
    sum = 0
    for elem in new_list:
        sum += elem
    moyenne = sum / len(new_list)
    return(moyenne)



def premier_pole_k(N, beta, angle_TT, angle_SDF, sens, sens_numerotation):
    """
    ATTENTION : Tous les angles sont à rentrer en degrés

    renvoie le premier pôle qui passe sur le capteur de flux en fonction
    des différents angles passés en paramètres
    """
    first = None
    val = None
    if sens_numerotation == 0:
        """
        sens trigo
        """
        rotation = [N-i for i in range(0, N)]
    if sens_numerotation == 1:
        """
        sens horaire
        """
        angle_SDF = 360-angle_SDF
        rotation = [i for i in range(1, N+1)]
    if beta < 0:
        beta += 360
    ecart_angulaire = 360 / N
    k = floor(beta / ecart_angulaire)
    X = beta - (k * ecart_angulaire)
    # sonde entre pôle 1+k et 2+k
    E = angle_SDF - angle_TT
    E_E = abs(E)
    if sens == 0:
        if E > 0 and E == ecart_angulaire:
            val = 2 + k
        elif E > 0 and E < ecart_angulaire:
            Y = ecart_angulaire - X
            if Y <= E:
                val = 2 + k
            elif Y > E:
                val = 1 + k
        elif E > 0 and E > ecart_angulaire:
            Y = ecart_angulaire - X
            l = floor(E / ecart_angulaire)
            W = E - (l * ecart_angulaire)
            if Y <= W:
                val = 2 + k + l
            elif Y > W:
                val = 1 + k + l
        elif E == 0:
            val = 1 + k
        elif E < 0 and E_E == ecart_angulaire:
            val = k
        elif E < 0 and E_E < ecart_angulaire:
            if X < E_E:
                val = k
            elif X >= E_E:
                val = 1 + k
        elif E < 0 and E_E > ecart_angulaire:
            l = floor(E_E / ecart_angulaire)
            W = E_E - (l * ecart_angulaire)
            if X < W:
                val = k - l
            elif X >= W:
                val = 1 + k - l
    elif sens == 1:
        if E > 0 and E == ecart_angulaire:
            if X == 0:
                val = 2 + k
            else:
                val = 3 + k
        elif E > 0 and E < ecart_angulaire:
            Y = ecart_angulaire - X
            if Y < E:
                val = 3 + k
            elif Y >= E:
                val = 2 + k
        elif E > 0 and E > ecart_angulaire:
            Y = ecart_angulaire - X
            l = floor(E / ecart_angulaire)
            W = E - (l * ecart_angulaire)
            if Y < W:
                val = 3 + k + l
            if Y >= W:
                val = 2 + k + l
        elif E == 0:
            if X == 0:
                val = 1 + k
            else:
                val = 2 + k
        elif E < 0 and E_E == ecart_angulaire:
            if X == 0:
                val = k
            else:
                val = k + 1
        elif E < 0 and E_E < ecart_angulaire:
            if X <= E_E:
                val = 1 + k
            elif X > E_E:
                val = 2 + k
        elif E < 0 and E_E > ecart_angulaire:
            l = floor(E_E / ecart_angulaire)
            W = E_E - (l * ecart_angulaire)
            if X <= W:
                val = 1 + k + l
            elif X > W:
                val = 2 + k + l
    """
    on peut obtenir une valeur > N ou négative donc on prend val et on applique modulo N
    """
    if val > N:
        first = val - N
    elif val < 1:
        first = val + N
    else:
        first = val
    real_first = 0
    if sens == 0:
        """
        sens trigo
        """
        real_first = first + 1
        if real_first > N:
            real_first -= N
    elif sens == 1:
        """
        sens horaire
        """
        real_first = first - 1
        if real_first < 1:
            real_first += N
    return real_first



def moyenne_par_pole(cretes, pole_actuel, sens_rotation, nb_poles):
    moyenne = [0 for _ in range(nb_poles)]
    values = [[] for _ in range(nb_poles)]
    if sens_rotation == 0:
        """
        sens trigo
        """
        for i in range(len(cretes)):
            values[pole_actuel - 1].append(cretes[i])
            pole_actuel -= 1
            if pole_actuel == 0:
                pole_actuel = nb_poles
    if sens_rotation == 1:
        """
        sens horaire
        """
        for i in range(len(cretes)):
            values[pole_actuel - 1].append(cretes[i])
            pole_actuel += 1
            if pole_actuel == nb_poles + 1:
                pole_actuel = 1
    for k in range(len(values)):
        sum = 0
        for l in range(len(values[k])):
            sum += values[k][l]
        moy = sum / len(values[k])
        moyenne[k] = moy
    return(moyenne)



def calcul_point_moyen(stockage):
    sum = 0
    nb_points = len(stockage)*len(stockage[0])
    for i in range(len(stockage)):
        for j in range(len(stockage[0])):
            sum += stockage[i][j]
    pt_moyen = sum / nb_points
    return(pt_moyen)



def normalisation_poles(moyenne_des_poles, pt_moyen):
    normalised = []
    for k in range(len(moyenne_des_poles)):
        normalised.append(moyenne_des_poles[k] / pt_moyen)
    return(normalised)



def main():
    """
    execute toutes les fonctions définies précédemment
    """
    print("Que souhaitez-vous étudier ? (pour fem tapez 0, pour entrefer tapez 1, et pour toptour tapez 2)")
    commande = int(input())
    nb_colonnes, nb_lignes, id_temps, les_fems, id_TT, les_entrefers = lecture_csv()
    del les_fems[-1]
    if commande == 1:
        """
        on s'occupe des entrefers
        """
        les_t = analyse_temps(nb_lignes, id_temps) # on lit la colonne du temps
        toptour = analyse_TT(nb_lignes, id_TT) # on lit la colonne du toptour
        seuil, sens = calcul_seuil(toptour, les_t) # on calcule le seuil du tt et le sens pour savoir comment on détecte
        les_tt = tour_machine(toptour, les_t, seuil, sens) # on renvoie les id de chaque tour machine(741, 1450, etc...)
        ecart_entre_t = les_t[1] - les_t[0]
        ecart_entre_tt = les_tt[1] - les_tt[0]
        tour_machine2 = ecart_entre_tt * ecart_entre_t
        paire_de_poles = round(tour_machine2 / 20)
        N = 2 * paire_de_poles # N = nombre de poles
        print("Veuillez entrer alpha_TT: l'angle du stator/de la partie fixe du toptour (en degrés)")
        angle_TT = round(float(input()))
        print("Veuillez entrer beta: l'angle du rotor/de la partie mobile du toptour (en degrés)")
        beta = round(float(input()))
        print("Dans quel sens tourne le rotor ? (entrez 0 pour sens trigonométrique et 1 pour sens horaire)")
        sens = int(input())
        print("Dans quel sens sont rentrés les angles des capteurs ? (entrez 0 pour sens trigonométrique et 1 pour sens horaire)")
        sens_numerotation = int(input())
        stockage = []
        """
        # truc
        entrefer_k = analyse_entrefer_k(nb_lignes, les_entrefers, les_entrefers[0])
        seuil2, sens2 = calcul_seuil(entrefer_k, les_t)
        print(seuil2)
        print(sens2)
        print(N)
        kk = selection_entrefer_k(entrefer_k, les_t, les_tt, seuil2, 1, N)
        print(kk)
        """
        for k in les_entrefers:
            entrefer_k = analyse_entrefer_k(nb_lignes, les_entrefers, k)
            seuil2, sens2 = calcul_seuil(entrefer_k, les_t)
            decoupage_entrefer_k = selection_entrefer_k(entrefer_k, les_t, les_tt, seuil2, 1, N)
            les_moyennes = []
            for t in decoupage_entrefer_k:
                moyenne =  calcul_valeur_moyenne(t, ecart_entre_t)
                les_moyennes.append(moyenne)
            print("Veuillez entrer angle_SDF_{}: l'angle de la sonde de flux {} (en degrés)".format(k-les_entrefers[0]+1, k-les_entrefers[0]+1))
            angle_SDF_k = int(input())
            first_pole_k = premier_pole_k(N, beta, angle_TT, angle_SDF_k, sens, sens_numerotation)
            print(first_pole_k)
            moyenne_des_poles = moyenne_par_pole(les_moyennes, first_pole_k, sens, N)
            print(moyenne_des_poles)
            stockage.append(moyenne_des_poles)
            """
            plt.subplot(321 + (k/2))
            plt.plot(les_t, entrefer_k)
            plt.title("les entrefers du capteur {}".format(k-6))
            plt.ylabel("les entrefers (mm)")
            plt.xlabel("temps (ms)")
            """
        les_numeros = [(l + 1) for l in range(len(stockage[0]))]
        pt_moyen = calcul_point_moyen(stockage)
        for i in range(len(stockage)):
            normalised = normalisation_poles(stockage[i], pt_moyen)
            plt.subplot(321 + i)
            plt.plot(les_numeros, normalised)
            plt.title("moyennes pour la sonde {}".format(i+1))
            plt.ylabel("les moyennes (V)")
            plt.xlabel("numéro de pôle")
        print(pt_moyen)
        


main()

plt.show()
