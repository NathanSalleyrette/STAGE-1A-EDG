# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 14:28:07 2020

@author: Nathan Salleyrette, modifié par André PETIT
"""

import csv
import pandas as pd

def lecture_csv(fichier):
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
    file = open(fichier,"r")
    # contenuCSV = io.TextIOBase.file.read()    # on lit la totalité du fichier sous forme de grande chaine de caractères
    test = csv.reader(file)
    for row in test:
        nb_colonnes = (len(row))
        nb_lignes += 1
    file.seek(0)   # revenir au début du fichier
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

def entetes_csv(fichier, colonnes):
    """
    retourne une liste avec le nom des entêtes du fichier CSV
    on fournit une liste des colonnes souhaitées
    """
    data = pd.read_csv(fichier, usecols=colonnes, encoding='ANSI')  # par défaut Panda lit de l'UTF-8
    return data.columns.values


def analyse_temps(fichier, nb_lignes, id_temps):
    """
    lit la colonne du temps, replace le point de depart à 0s et convertit en ms
    pour plus de lisibilité et renvoie une liste avec tous les instants t
    """

    t = []
    file=open(fichier,"r")
    ligne = file.readline()
 #   file = io.StringIO(contenuCSV)  # on ne relit pas le fichier, mais son image
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

def analyse_TT(fichier, nb_lignes, id_TT):
    """
    lit la colonne qui contient des infos sur le toptour et replace toutes
    ces valeurs dans une liste toptour
    """

    toptour = []
    file=open(fichier,"r")
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


def analyse_fem_k(fichier, nb_lignes, les_fems, k):
    """
    lit la colonne qui contient des infos sur les fems et replace toutes
    les valeurs dans une liste fem_k
    """
    toptour = []
    fem_k = []
    file=open(fichier,"r")
    ligne = file.readline()
    # on lit la 1ere ligne
    # on lit les autres lignes de données
    for i in range(2, nb_lignes + 1):
        ligne = file.readline()
        # pour la fem k
        element2 = ligne.split(',')[k]
        fe2 = float(element2)
        fem_k.append(fe2)
    file.close()
    return fem_k


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
