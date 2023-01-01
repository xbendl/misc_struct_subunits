#!/usr/bin/python
# -*- coding: utf-8 -*-

# PBI - Projekt, varianta B
#
# strElement.py - Modul s definici tridy pro ulozeni informace o strukturnim
#                 elementu shodnem pro vsechny vysetrovane proteinove struktury
#
# Autor: Jaroslav Bendl, xbendl00
#

import sys

class StrElement:
    """Trida uchovavajici informace o strukturnim elementu."""

    def __init__(self, argAminoAcidPair1, argAminoAcidPair2, positionList):
        self.aminoAcidPair1 = argAminoAcidPair1    # Prvni par sousedicich aminokyselin na kostre proteinu
        self.aminoAcidPair2 = argAminoAcidPair2    # Druhy par sousedicich aminokyselin na kostre proteinu
        self.proteinPosition = []                   # List pozic paru aminokyselin pro konkretni proteinove struktury

        if(len(positionList) > 0):
            self.proteinPosition = positionList    
        

    def addProteinPosition(self, proteinName, i1, i2, j1, j2):
        """Vklada pozici vyskytu strukturniho elementu pro danou proteinovou strukturu."""
        self.proteinPosition.append([proteinName, i1, i2, j1, j2])
        

    def getStructureElement(self):
        """Vrati nalezeny strukturni element vcetne pozic vyskytu pro jednotlive proteinove struktury."""
        
        elementReport = self.aminoAcidPair1 + " " + self.aminoAcidPair2 + " "
        
        for i in range(len(self.proteinPosition)):
            elementReport += self.proteinPosition[i][0] + ":" + str(self.proteinPosition[i][1]) + ":" + str(self.proteinPosition[i][2]) + ":" + str(self.proteinPosition[i][3]) + ":" + str(self.proteinPosition[i][4]) + " "
        
        return elementReport
        
  