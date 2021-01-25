#!/usr/bin/env python3
import csv
import numpy as np
from scipy import integrate as it
from scipy.integrate import simps
from scipy import signal
from math import *
import matplotlib.pyplot as plt
import dbManip
import traitement



def main():
    """
    execute toutes les fonctions définies précédemment
    """
    fichier =  "Auzat_new.csv"
    print("Que souhaitez-vous étudier ? (pour fem tapez 0, pour entrefer tapez 1, et pour toptour tapez 2)")
    commande = int(input())
    print("Est-ce-qu'on normalise les résultats ? (pour oui tapez 0, pour non tapez 1)")
    normalize = int(input())
    nb_colonnes, nb_lignes, id_temps, les_fems, id_TT, les_entrefers = dbManip.lecture_csv(fichier)
    del les_fems[-1]
    if commande == 0:
        """
        on s'occupe des fem
        """
        les_t = dbManip.analyse_temps(fichier, nb_lignes, id_temps) # on lit la colonne du temps
        toptour = dbManip.analyse_TT(fichier, nb_lignes, id_TT) # on lit la colonne du toptour
        seuil, sens = traitement.calcul_seuil(toptour, les_t) # on calcule le seuil du tt et le sens pour savoir comment on détecte
        les_tt = traitement.tour_machine(toptour, les_t, seuil, sens) # on renvoie les id de chaque tour machine(741, 1450, etc...)
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
        for k in les_fems:
            #on chope temps, fem_k, toptour
            fem_k = dbManip.analyse_fem_k(fichier, nb_lignes, les_fems, k)
            new_fem_k, new_t = traitement.retirer_offset(fem_k, les_t, les_tt, tour_machine2)
            new_fem_int_k = it.cumtrapz(new_fem_k, new_t, initial=0)
            last_fem_k, last_t = traitement.recalculer_offset(new_fem_int_k, new_t, les_tt, tour_machine2)
            t2_k = [0 for _ in range(len(last_fem_k))]
            for i in range(len(last_fem_k)):
                t2_k[i] = abs(last_fem_k[i]) # on redresse la tension
            les_cretes_k = traitement.detection_cretes(last_t, t2_k, len(last_t))
            print("Veuillez entrer angle_SDF_{}: l'angle de la sonde de flux {} (en degrés)".format(k, k))
            angle_SDF_k = round(float(input()))
            first_pole_k = traitement.premier_pole_k(N, beta, angle_TT, angle_SDF_k, sens, sens_numerotation)
            print("Le premier pôle détecté est le pôle {}".format(first_pole_k))
            moyenne_des_poles = traitement.moyenne_par_pole(les_cretes_k, first_pole_k, sens, N)
            print(moyenne_des_poles)
            stockage.append(moyenne_des_poles)
        les_numeros = [(l + 1) for l in range(len(stockage[0]))]
        if normalize == 0:
            pt_moyen = traitement.calcul_point_moyen(stockage)
            for i in range(len(stockage)):
                normalised = traitement.normalisation_poles(stockage[i], pt_moyen)
                plt.subplot(322 + i)
                plt.plot(les_numeros, normalised)
                plt.title("moyennes pour la sonde {}".format(i + 1))
                plt.ylabel("valeurs relatives")
                plt.xlabel("numéro de pôle")
            # cacher le bloc a venir si on veut les 4 courbes différentes
            plt.subplot(321)
            for j in range(len(stockage)):
                plt.plot(les_numeros, stockage[j], label = 'spire {}'.format(j + 1))
                plt.title("les courbes superposées")
                plt.ylabel("les moyennes (Wb)")
                plt.xlabel("numéro de pôle")
                plt.legend()
            print("le flux moyen vaut {} Wb".format(pt_moyen))
        elif normalize == 1:
            for i in range(len(stockage)):
                plt.subplot(322 + i)
                plt.plot(les_numeros, stockage[i])
                plt.title("moyenne pour la sonde {}".format(i + 1))
                plt.ylabel("les moyennes (Wb)")
                plt.xlabel("numéro de pôle")
        
        data = []
        data1 = ['X']
        data2 = ['spire {}'.format(l+1) for l in range(len(stockage))]
        data1.extend(data2)
        data.append(data1)
        for k in range(N):
            data3 = ['pole {}'.format(k+1)]
            data4 = [str(stockage[m][k]) for m in range(len(stockage))]
            data3.extend(data4)
            data.append(data3)

        with open('data.csv', mode='w') as file:
            writer = csv.writer(file, delimiter=';')
            # write data
            for ligne in data:
                writer.writerow(ligne)
    if commande == 1:
        """
        on s'occupe des entrefers
        """
        les_t = dbManip.analyse_temps(fichier, nb_lignes, id_temps) # on lit la colonne du temps
        toptour = dbManip.analyse_TT(fichier, nb_lignes, id_TT) # on lit la colonne du toptour
        seuil, sens = traitement.calcul_seuil(toptour, les_t) # on calcule le seuil du tt et le sens pour savoir comment on détecte
        les_tt = traitement.tour_machine(toptour, les_t, seuil, sens) # on renvoie les id de chaque tour machine(741, 1450, etc...)
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
        for k in les_entrefers:
            entrefer_k = dbManip.analyse_entrefer_k(nb_lignes, les_entrefers, k)
            seuil2, sens2 = traitement.calcul_seuil(entrefer_k, les_t)
            decoupage_entrefer_k = traitement.selection_entrefer_k(entrefer_k, les_t, les_tt, seuil2, 1, N)
            les_moyennes = []
            for t in decoupage_entrefer_k:
                moyenne =  traitement.calcul_valeur_moyenne(t, ecart_entre_t)
                les_moyennes.append(moyenne)
            print("Veuillez entrer angle_capteur_{}: l'angle du capteur entrefer n°{} (en degrés)".format(k-les_entrefers[0]+1, k-les_entrefers[0]+1))
            angle_SDF_k = int(input())
            first_pole_k = traitement.premier_pole_k(N, beta, angle_TT, angle_SDF_k, sens, sens_numerotation)
            print("Le premier pôle détecté est le pôle {}$".format(first_pole_k))
            moyenne_des_poles = traitement.moyenne_par_pole(les_moyennes, first_pole_k, sens, N)
            print(moyenne_des_poles)
            stockage.append(moyenne_des_poles)
        les_numeros = [(l + 1) for l in range(len(stockage[0]))]
        if normalize == 0:
            pt_moyen = traitement.calcul_point_moyen(stockage)
            for i in range(len(stockage)):
                normalised = traitement.normalisation_poles(stockage[i], pt_moyen)
                plt.subplot(322 + i)
                plt.plot(les_numeros, normalised)
                plt.title("moyennes pour la sonde {}".format(i+1))
                plt.ylabel("les moyennes relatives")
                plt.xlabel("numéro de pôle")
            # cacher le bloc a venir si on veut les 4 courbes différentes
            plt.subplot(321)
            for j in range(len(stockage)):
                plt.plot(les_numeros, stockage[j], label = 'spire {}'.format(j + 1))
                plt.title("les courbes superposées")
                plt.ylabel("les moyennes (mm)")
                plt.xlabel("numéro de pôle")
                plt.legend()
            print("L'entrefer moyen vaut {} mm".format(pt_moyen))
        elif normalize == 1:
            for i in range(len(stockage)):
                plt.subplot(322 + i)
                plt.plot(les_numeros, stockage[i])
                plt.title("moyenne pour la sonde {}".format(i + 1))
                plt.ylabel("les moyennes (mm)")
                plt.xlabel("numéro de pôle")
                
        data = []
        data1 = ['X']
        data2 = ['capteur entrefer n°{}'.format(l+1) for l in range(len(stockage))]
        data1.extend(data2)
        data.append(data1)
        for k in range(N):
            data3 = ['pole {}'.format(k+1)]
            data4 = [str(stockage[m][k]) for m in range(len(stockage))]
            data3.extend(data4)
            data.append(data3)

        with open('data.csv', mode='w') as file:
            writer = csv.writer(file, delimiter=';')
            # write data
            for ligne in data:
                writer.writerow(ligne)
    if commande == 2:
        """
        on s'occupe du toptour
        """
        les_t = dbManip.analyse_temps(fichier, nb_lignes, id_temps) # on lit la colonne du temps
        toptour = dbManip.analyse_TT(fichier, nb_lignes, id_TT) # on lit la colonne du toptour
        plt.plot(les_t, toptour)
        plt.title("toptour au cours du temps")
        plt.ylabel("toptour (en micrometres)")
        plt.xlabel("temps (ms)")


main()

plt.show()
