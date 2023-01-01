#!/usr/bin/python
# -*- coding: utf-8 -*-

# PBI - Projekt, varianta B
#
# strUnit.py - Modul s definici tridy pro ulozeni informace o nalezenem 
#             strukturnim celku.
#
# Autor: Jaroslav Bendl, xbendl00
#

import sys
from strElement import StrElement

class StrUnit:
    """Trida uchovavajici informace o nalezenem stukturnim celku."""

    def __init__(self):
        self.aminoAcidSeq = ""   # vypocitana sekvence aminokyselin
        self.aminoAcidPos = []   # vypocitane pozice aminokyselin v jednotlivych proteinech
        self.element = []        # seznam strukturnich jednotek
        self.elementIndex = []   # seznam indexu pouzitych strukturnich jednotek z puvodniho seznamu
        self.sealed = False      # priznak urcujici, zda-li ma jeste smysl provadet pokusy o rozsirovani strukturniho celku
        
        
    def addElement(self, argElement, argElementIndex, newAminoAcidIndex):
        """Pokusi se pridat element do seznamu strukturnich elementu (overuje topologicke usporadani v proteinech)."""
        #print "N: " + str(newAminoAcidIndex)
        # Pokud je seznam strukturnich elementu prazdny, prida se vzdy
        if(len(self.element) == 0):
            self.aminoAcidSeq = argElement.aminoAcidPair1 + argElement.aminoAcidPair2
            for i in range(len(argElement.proteinPosition)):
                subList = argElement.proteinPosition[i][:]
                self.aminoAcidPos.append(subList)
            
            self.element.append(argElement)
            self.elementIndex.append(argElementIndex)
            return True
        
        # Overeni, zda-li jiz nahodou nejsou vsechny pozice v ramci daneho celku zarazeny
        includeCount = 0
        isIncluded = False
        for j in range(len(self.aminoAcidPos[0]) - 1):
            if((self.aminoAcidPos[0][j+1] == argElement.proteinPosition[0][1])  or \
                (self.aminoAcidPos[0][j+1] == argElement.proteinPosition[0][2]) or \
                (self.aminoAcidPos[0][j+1] == argElement.proteinPosition[0][3]) or \
                (self.aminoAcidPos[0][j+1] == argElement.proteinPosition[0][4])):
                    includeCount = includeCount + 1
        if(includeCount == 4):
            isIncluded = True

        # Pokud neni seznam strukturnich elementu prazdny, je nutne overit topologicke usporadani v proteinech
        index = 1       # index pro vlozeni "nove" aminokyseliny (topologicka spravnost = stejny pro kazdou strukturu)
        fixIndex = -1   # pomocna - pro kontrolu, zda se nemeni index "nove" aminokyseliny
        
        # Kontrola pres vsechny struktury
        for i in range(len(self.aminoAcidPos)):

            # Hledani indexu pro zarazeni "nove" aminokyseliny
            index = 0
            #print len(self.aminoAcidPos[i])
            #print newAminoAcidIndex
            #print str(len(argElement.proteinPosition[i])) + " - " + str(len(self.aminoAcidPos[i]))
            while index < (len(self.aminoAcidPos[i]) - 1):
                #print str(argElement.proteinPosition[i][newAminoAcidIndex]) + " > " + str(self.aminoAcidPos[i][index+1])
                #print str(argElement.proteinPosition[0][newAminoAcidIndex]) + " > " + str(self.aminoAcidPos[0][1])
                if(argElement.proteinPosition[i][newAminoAcidIndex] > self.aminoAcidPos[i][index+1]):
                    index = index + 1
                else:
                    break
            index = index + 1
            if(fixIndex == -1):    # pocatecni nastaveni pomocne promenne fixIndex
                fixIndex = index
            
            if(fixIndex != index): # odlisne indexy pro zarazeni "nove" aminokyseliny => neni topologicky spravne => konec
                return False
            
            if(isIncluded == False):
                self.aminoAcidPos[i].insert(index, argElement.proteinPosition[i][newAminoAcidIndex])
        
        # Pridani "nove" aminokyseliny do textove reprezentace
        if(isIncluded == False):
            if(newAminoAcidIndex <= 2):
                self.aminoAcidSeq = self.aminoAcidSeq[0:(index-1)] + argElement.aminoAcidPair1[newAminoAcidIndex - 1] + self.aminoAcidSeq[(index-1):]
            else:
                self.aminoAcidSeq = self.aminoAcidSeq[0:(index-1)] + argElement.aminoAcidPair2[newAminoAcidIndex - 3] + self.aminoAcidSeq[(index-1):]
            
        # Pridani nove strukturni jednotky do strukturniho celku a zapsani jeho indexu
        self.element.append(argElement)
        i = 0
        while i < len(self.elementIndex):
            if(argElementIndex > self.elementIndex[i]):
                i = i + 1
            else:
                break
        self.elementIndex.insert(i, argElementIndex)
        
        return True
