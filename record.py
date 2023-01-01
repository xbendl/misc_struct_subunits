#!/usr/bin/python
# -*- coding: utf-8 -*-

# PBI - Projekt, varianta B
#
# record.py - Modul s definici tridy pro ulozeni informaci o vysetrovanych
#             proteinech a jejich umisteni na disku
#
# Autor: Jaroslav Bendl, xbendl00
#

class Record:
    """Trida uchovavajici informace o vysetrovanych proteinech a jejich umisteni na disku"""
    
    def __init__(self, argProteinId, argProteinFile, argProteinChain):
        self.proteinId = argProteinId         # Identifikator proteinu v databazi struktur PDB
        self.proteinFile  = argProteinFile    # Cesta k zaznamu struktury proteinu (ziskane z databaze PDB)
        self.proteinChain = argProteinChain   # Retezec na strukture, ktery nas zajima (empty => struktura ma pouze 1 retezec)
