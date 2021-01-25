# -*- coding: utf-8 -*-
from scipy.integrate import simps
from math import *




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
        seuil = id_max + (ecart_max * 2 / 3)
        sens = 1
    elif len(dessus) < len(dessous):
        seuil = id_max - (ecart_max * 2 / 3)
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



def retirer_offset(fem_k, les_t, les_tt, tour_machine2):
    """
    calcule l'aire entre chaque tour machine et retranche cet offset au signal
    il renvoie une liste de fem sans offset et une liste de temps associé
    à cette même fem
    """
    new_fem = []
    new_t = []
    for j in range(len(les_tt) - 1):
        offset = (calcul_aire(les_t, fem_k, les_tt[j], les_tt[j + 1])) / tour_machine2
        for m in range(les_tt[j], les_tt[j + 1]):
            new_y = (fem_k[m] - offset)
            new_fem.append(new_y)
            new_t.append(les_t[m])
    new_fem.append(fem_k[les_tt[-1]] - offset)
    new_t.append(les_t[les_tt[-1]])
    return(new_fem, new_t)


def calcul_aire(temps, fem, depart, arrivee):
    """
    calcule l'aire entre un indice de depart et un indice d'arrivée
    """
    les_y = [0 for _ in range(arrivee - depart + 1)]
    les_x = [0 for _ in range(arrivee - depart + 1)]
    for i in range(arrivee - depart + 1):
        les_x[i] = temps[depart + i]
        les_y[i] = fem[depart + i]
    area = simps(les_y, les_x)
    return area



def recalculer_offset(new_fem_int_k, new_t, les_tt, tour_machine2):
    """
    calcule à nouveau l'aire entre chaque tour machine et retranche cet offset
    au signal intégré. Il renvoie une liste de fem intégrée sans offset et
    la liste de temps associée à ce signal
    """
    last_fem = []
    last_t = []
    first_tt = les_tt[0]
    periodes = [(l - first_tt) for l in les_tt]
    for j in range(len(les_tt) - 1):
        offset = (calcul_aire(new_t, new_fem_int_k, periodes[j], periodes[j + 1])) / tour_machine2
        for m in range(periodes[j], periodes[j + 1]):
            last_y = (new_fem_int_k[m] - offset)
            last_fem.append(last_y)
            last_t.append(new_t[m])
    return(last_fem, last_t)


def detection_cretes(t, fem, nb_points):
    """
    renvoie les maximums de chaque crêtes composant un signal de fem en
    fonction du temps qui contient nb_points
    """
    cretes = []
    clef = 1
    seuil = 0
    max = 0
    for i in range(nb_points):
        if fem[i] > max and clef == 1:
            max = fem[i]
            seuil = max / 4
        if fem[i] < seuil and clef == 1:
            cretes.append(max)
            clef = 0
        if fem[i] > seuil and clef == 0:
            max = fem[i]
            clef = 1
    if cretes[0] < (0.9*cretes[1]):
        del cretes[0] # on supprime la 1ere valeur si elle est absurde
    return cretes



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
    """
    renvoie la valeur moyenne des valeurs des différentes crêtes
    associées au pole_actuel
    """
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
    """
    calcul le point moyen d'une matrice nommée stockage
    """
    sum = 0
    nb_points = len(stockage)*len(stockage[0])
    for i in range(len(stockage)):
        for j in range(len(stockage[0])):
            sum += stockage[i][j]
    pt_moyen = sum / nb_points
    return(pt_moyen)



def normalisation_poles(moyenne_des_poles, pt_moyen):
    """
    normalise les points de la liste moyenne_des_poles en divisant
    tous ces points par pt_moyen
    """
    normalised = []
    for k in range(len(moyenne_des_poles)):
        normalised.append(moyenne_des_poles[k] / pt_moyen)
    return(normalised)


def selection_entrefer_k(entrefer_k, les_t, les_tt, seuil, sens, N):
    """
    sélectionne les points au sommet de chaque crête dans le signal entrefer 
    en fonction du temps
    """
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
    """
    calcule la valeur moyenne au sommet de chaque crête du signal
    entrefer en fonction du temps
    """
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
