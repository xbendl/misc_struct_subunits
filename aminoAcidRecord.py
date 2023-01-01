#!/usr/bin/python
# -*- coding: utf-8 -*-

# PBI - Projekt, varianta B
#
# aminoAcidReport.py - Modul s definici tridy pro ulozeni informace o
#                      aminokyseline na pateri pozadovaneho rezetce vysetrovane
#                      proteinove struktury z PDB databaze.
#
# Autor: Jaroslav Bendl, xbendl00
#

import math

# Tabulka pro konverzi mezi tripismenou a jednopismenou zkratkou aminokyseliny
AMINO = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 
         'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
         'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
         'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
         'SEC': 'U', 'PYL': 'O',}  

class AminoAcidRecord:
    """Trida uchovavajici informace o aminokyseline na pateri pozadovaneho retezce vysetrovane struktury z PDB databaze."""

    def __init__(self, argResiNumber, argAminoAcid, argCoordX, argCoordY, argCoordZ, angstromDist):
        self.resiNumber = argResiNumber    # Poradi aminokyseliny na kostre
        self.aminoAcid  = self.__convertAminoAcidFormat(argAminoAcid)     # Reprezentovana aminokyselina
        self.coordX = argCoordX            # X-souradnice C-alpha uhliku dane aminokyseliny
        self.coordY = argCoordY            # Y-souradnice C-alpha uhliku dane aminokyseliny
        self.coordZ = argCoordZ            # Z-souradnice C-alpha uhliku dane aminokyseliny
        self.neighbourDistance = float(-1)
        self.neighbourDistanceLimit = angstromDist
        
        
    def __convertAminoAcidFormat(self, argAminoAcid):
        """Konvertuje oznaceni aminokyseliny z tripismenoveho kodu do jednopismenoveho."""
        return AMINO[argAminoAcid]
    
    
    def setNeighbourDistance(self, neighbour):
        """Nastavi vzdalenost k sousedovi na pateri retezce proteinove struktury."""
        self.neighbourDistance = float(math.sqrt(pow((neighbour.coordX - self.coordX), 2) + pow((neighbour.coordY - self.coordY), 2) + pow((neighbour.coordZ - self.coordZ), 2)))


    def isNearNeighbour(self):
        """Zjisti zda vzdalenost aktualniho c-alpha uhliku od sousedni aminokyseliny je v povolenem limitu."""
        return (self.neighbourDistance < self.neighbourDistanceLimit)
        
        
    def isNear(self, neighbour):
        """Zjisti zda vzdalenost aktualniho c-alpha uhliku od urcene aminokyseliny je v povolenem limitu."""
        distance = float(math.sqrt(pow((neighbour.coordX - self.coordX), 2) + pow((neighbour.coordY - self.coordY), 2) + pow((neighbour.coordZ - self.coordZ), 2)))
        #print distance
        return (distance < self.neighbourDistanceLimit)
        
