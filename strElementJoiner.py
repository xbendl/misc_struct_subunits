#!/usr/bin/python
# -*- coding: utf-8 -*-

# PBI - Projekt, varianta B
#
# strElementJoiner.py - Modul pro identifikaci mnoziny nejvetsich vzajemne
#                       nekolidujicich strukturnich celku.
#
# Autor: Jaroslav Bendl, xbendl00
#

import copy
from strElement import StrElement
from strUnit import StrUnit
from pymol import cmd

class StrElementJoiner:
    """Ziska seznam strukturnich celku shodnych pro danou sadu proteinovych struktur."""

    def __init__(self):
        self.unitList = []  # seznam vsech strukturnich celku (jeste se nezkouma, zda-li koliduji ci nikoliv)

    def getStructureUnits(self, elementList):
        """Ziska seznam strukturnich celku shodnych pro danou sadu proteinovych struktur."""
       
        # Generovani strukturniho celku je vyzkouseno z kazdeho elementu
        for i in range(len(elementList)):
            unit = StrUnit()
            unit.addElement(elementList[i], i, -1)

            # Hledani navazujiciho elementu v seznamu vsech elementu
            for j in range(len(elementList)):
                if(j == i):
                    continue
                
                # Priblizne stanoveni, zda jsou elementy navazujici (podle nazvu aminokyselin - 3/4 musi byt shodne)
                cond1 = (unit.element[0].aminoAcidPair1 == elementList[j].aminoAcidPair1) and \
                  ((unit.element[0].aminoAcidPair2[0] == elementList[j].aminoAcidPair2[1]) or \
                  (unit.element[0].aminoAcidPair2[1] == elementList[j].aminoAcidPair2[0]))
                cond2 = (unit.element[0].aminoAcidPair2 == elementList[j].aminoAcidPair2) and \
                  ((unit.element[0].aminoAcidPair1[0] == elementList[j].aminoAcidPair1[1]) or \
                  (unit.element[0].aminoAcidPair1[1] == elementList[j].aminoAcidPair1[0]))

                # Presne stanoveni, zda jsou elementy navazujici (overeni pozic - 3/4 musi byt shodne)
                if(cond1 or cond2):
                    c1 = c2 = c3 = c4 = 0
                    for k in range(len(unit.element[0].proteinPosition)):
                        c1 = c2 = c3 = c4 = 0
                        
                        # Detekce vztahu elementu typu AB CD => XA CD
                        if((unit.element[0].proteinPosition[k][3]  == elementList[j].proteinPosition[k][3]) and \
                            (unit.element[0].proteinPosition[k][4] == elementList[j].proteinPosition[k][4]) and \
                            (unit.element[0].proteinPosition[k][1] == elementList[j].proteinPosition[k][2]) and \
                            (unit.element[0].proteinPosition[k][2] != elementList[j].proteinPosition[k][1])):
                            c1 = 1
                            continue
                        else:
                            if(c1 == 1):
                                c1 = c2 = c3 = c4 = 0
                                break
                        # Detekce vztahu elementu typu AB CD => BX CD
                        if((unit.element[0].proteinPosition[k][3]  == elementList[j].proteinPosition[k][3]) and \
                            (unit.element[0].proteinPosition[k][4] == elementList[j].proteinPosition[k][4]) and \
                            (unit.element[0].proteinPosition[k][1] != elementList[j].proteinPosition[k][2]) and \
                            (unit.element[0].proteinPosition[k][2] == elementList[j].proteinPosition[k][1])):
                            c2 = 2
                            continue
                        else:
                            if(c2 == 2):
                                c1 = c2 = c3 = c4 = 0
                                break
                        # Detekce vztahu elementu typu AB CD => AB XC
                        if((unit.element[0].proteinPosition[k][1]  == elementList[j].proteinPosition[k][1]) and \
                            (unit.element[0].proteinPosition[k][2] == elementList[j].proteinPosition[k][2]) and \
                            (unit.element[0].proteinPosition[k][3] == elementList[j].proteinPosition[k][4]) and \
                            (unit.element[0].proteinPosition[k][4] != elementList[j].proteinPosition[k][3])):
                            c3 = 3
                            continue
                        else:
                            if(c3 == 3):
                                c1 = c2 = c3 = c4 = 0
                                break
                        # Detekce vztahu elementu typu AB CD => AB DX
                        if((unit.element[0].proteinPosition[k][1]  == elementList[j].proteinPosition[k][1]) and \
                            (unit.element[0].proteinPosition[k][2] == elementList[j].proteinPosition[k][2]) and \
                            (unit.element[0].proteinPosition[k][3] != elementList[j].proteinPosition[k][4]) and \
                            (unit.element[0].proteinPosition[k][4] == elementList[j].proteinPosition[k][3])):
                            c4 = 4
                            continue
                        else:
                            if(c4 == 4):
                                c1 = c2 = c3 = c4 = 0
                                break
                    
                    index = c1 + c2 + c3 + c4  # urcuje index predstavujici "novou" aminokyselinu v strukturni jednotce
                    
                    # Strukturni element je navazujici (v jednom ze 4 moznosti)
                    if(index > 0):                       
                        # Kontrola topologie
                        if(unit.addElement(elementList[j], j, index) == True):
                            # Kontrola, zda jiz neexistuje struktura se zcela totoznymi strukturnimi elementy
                            breakFlag = False
                            for k in range(len(self.unitList)):
                                if(self.unitList[k].elementIndex == unit.elementIndex):
                                    breakFlag = True
                            if(breakFlag == False):
                                self.unitList.append(unit)
                                #print "Strukturni jednotky: " + unit.aminoAcidSeq

                        # Reset pomocne struktury pro novy strukturni celek
                        unit = StrUnit()
                        unit.addElement(elementList[i], i, -1)

        # Opakovane rozsirovani strukturnich celku az do rozgenerovani vsech variant 
        changeFlag = True
        while(changeFlag == True):
            i = 0
            changeFlag = False
            
            # Rozsirovani strukturnich celku o jeden strukturni element
            while(i < len(self.unitList)):
                # Rozsirovani probiha pouze u dosud nerozsirovanych strukturnich jednotek
                if(self.unitList[i].sealed == False):
                    # Pokus o rozsireni se zjistenim, zda-li k rozsireni doslo ci nikoliv
                    if(self.extendStructureUnit(elementList, self.unitList[i]) == True):
                        changeFlag = True

                i = i + 1
        
        # Kompletni vypis
        #for i in range(len(self.unitList)):
        #    print self.unitList[i].aminoAcidSeq
            
                
    def extendStructureUnit(self, elementList, unit):
        """Pokusi se rozsirit dany strukturni celek o dalsi strukturni element."""
        changeFlag = False
        #raw_input()
        
        # Hledani navazujiciho elementu v seznamu vsech elementu
        for i in range(len(elementList)):
            # Preskoceni prvku, ktere uz ve strukturnim elementu jsou
            continueFlag = False
            for j in range(len(unit.elementIndex)): 
                if(unit.elementIndex[j] == i):
                    continueFlag = True
                    break
            if(continueFlag == True):
                continue
                        
            # Zkouseni pripojeni pro kazdy strukturni element v ramci strukturniho celku
            j = 0
            for j in range(len(unit.element)):
                # Priblizne stanoveni, zda je element navazujici
                cond1 = (unit.element[j].aminoAcidPair1 == elementList[i].aminoAcidPair1)    and \
                  ((unit.element[j].aminoAcidPair2[0]   == elementList[i].aminoAcidPair2[1]) or \
                  (unit.element[j].aminoAcidPair2[1]    == elementList[i].aminoAcidPair2[0]))
                cond2 = (unit.element[j].aminoAcidPair2 == elementList[i].aminoAcidPair2)    and \
                  ((unit.element[j].aminoAcidPair1[0]   == elementList[i].aminoAcidPair1[1]) or \
                  (unit.element[j].aminoAcidPair1[1]    == elementList[i].aminoAcidPair1[0]))
            
                if((cond1 == True) or (cond2 == True)):
                    break

            # Presne stanoveni, zda jsou elementy navazujici (overeni pozic - 3/4 musi byt shodne)
            if(cond1 or cond2):
                c1 = c2 = c3 = c4 = 0
                for k in range(len(unit.element[0].proteinPosition)):
                    c1 = c2 = c3 = c4 = 0
                    
                    # Detekce vztahu elementu typu AB CD => XA CD
                    if((unit.element[j].proteinPosition[k][3]  == elementList[i].proteinPosition[k][3]) and \
                        (unit.element[j].proteinPosition[k][4] == elementList[i].proteinPosition[k][4]) and \
                        (unit.element[j].proteinPosition[k][1] == elementList[i].proteinPosition[k][2]) and \
                        (unit.element[j].proteinPosition[k][2] != elementList[i].proteinPosition[k][1])):
                        c1 = 1
                        continue
                    else:
                        if(c1 == 1):
                            c1 = c2 = c3 = c4 = 0
                            break
                    # Detekce vztahu elementu typu AB CD => BX CD
                    if((unit.element[j].proteinPosition[k][3]  == elementList[i].proteinPosition[k][3]) and \
                        (unit.element[j].proteinPosition[k][4] == elementList[i].proteinPosition[k][4]) and \
                        (unit.element[j].proteinPosition[k][1] != elementList[i].proteinPosition[k][2]) and \
                        (unit.element[j].proteinPosition[k][2] == elementList[i].proteinPosition[k][1])):
                        c2 = 2
                        continue
                    else:
                        if(c2 == 2):
                            c1 = c2 = c3 = c4 = 0
                            break
                    # Detekce vztahu elementu typu AB CD => AB XC
                    if((unit.element[j].proteinPosition[k][1]  == elementList[i].proteinPosition[k][1]) and \
                        (unit.element[j].proteinPosition[k][2] == elementList[i].proteinPosition[k][2]) and \
                        (unit.element[j].proteinPosition[k][3] == elementList[i].proteinPosition[k][4]) and \
                        (unit.element[j].proteinPosition[k][4] != elementList[i].proteinPosition[k][3])):
                        c3 = 3
                        continue
                    else:
                        if(c3 == 3):
                            c1 = c2 = c3 = c4 = 0
                            break
                    # Detekce vztahu elementu typu AB CD => AB DX
                    if((unit.element[j].proteinPosition[k][1]  == elementList[i].proteinPosition[k][1]) and \
                        (unit.element[j].proteinPosition[k][2] == elementList[i].proteinPosition[k][2]) and \
                        (unit.element[j].proteinPosition[k][3] != elementList[i].proteinPosition[k][4]) and \
                        (unit.element[j].proteinPosition[k][4] == elementList[i].proteinPosition[k][3])):
                        c4 = 4
                        continue
                    else:
                        if(c4 == 4):
                            c1 = c2 = c3 = c4 = 0
                            break
                
                index = c1 + c2 + c3 + c4  # urcuje index predstavujici "novou" aminokyselinu v strukturni jednotce

                # Strukturni element je navazujici (v jednom ze 4 moznosti)
                if(index > 0):
                    # Kontrola topologie
                    potentialUnit = copy.deepcopy(unit)
                    continueFlag = False
                    if(potentialUnit.addElement(elementList[i], i, index) == True):
                        for k in range(len(self.unitList)):
                            if(potentialUnit.aminoAcidPos == self.unitList[k].aminoAcidPos):
                                continueFlag = True
                        if(continueFlag == True):
                            continue
                        self.unitList.append(potentialUnit)
                        changeFlag = True
                        changeFlag = True
                        
        unit.sealed = True;  # priznak ze u daneho strukturniho elementu jiz nema probihat rozsirovani  
        
        return changeFlag        


    def getBestStructureUnits(self, unitCount):
        """Ulozi nalezene nejvetsi vzajemne nekolidujici strukturni celky do vystupniho souboru."""
        
        bestUnitList = []  # seznam nejvestich vzajemne nekolidujicich strukturnich celku
        
        # Hledani nejvetsiho nekolidujiciho strukturniho celku je opakovano do pozadovaneho poctu celku celku
        while((len(bestUnitList) < unitCount) and (len(self.unitList) > 0)):
            bestUnit = -1        # delka nejvetsiho strukturniho celku
            bestUnitIndex = -1   # index pro pristup k nejvetsimu strukturnimu celku

            # Prochazi se cely seznam strukturnich celku
            i = 0
            while i < len(self.unitList):
                # Pokud je aktualni strukturni celek vetsi, poznamena se jeho velikost a index
                if(len(self.unitList[i].aminoAcidSeq) > bestUnit):
                    bestUnit = len(self.unitList[i].aminoAcidSeq)
                    bestUnitIndex = i
                i = i + 1
            
            # Pokud je nalezeny strukturni celek nekolidujici s ostatnimi nejvetsimi, je k nim zarazen
            if(self.detectColision(bestUnitList, self.unitList[bestUnitIndex]) == False):
                bestUnitList.append(self.unitList[bestUnitIndex])
            
            # Nejvetsi strukturni celek je ze seznamu odstranen (at jiz koliduje ci nikoliv)
            self.unitList.pop(bestUnitIndex)
        
        return bestUnitList


    def detectColision(self, bestUnitList, unit):
        """Detekuje kolizi strukturniho celku s mnozinou strukturnich celku."""
        
        # Kontrola, zda vzajemne nekoliduji pozice nejvetsiho strukturniho celku s jakymkoliv jinym z mnoziny nejvetsich struktrunich celku
        for i in range(len(bestUnitList)):
            for j in range(len(bestUnitList[i].aminoAcidPos)):
                for k in range(len(bestUnitList[i].aminoAcidPos[j]) - 1):
                    for l in range(len(unit.aminoAcidPos[j]) - 1):
                        if(bestUnitList[i].aminoAcidPos[j][k+1] == unit.aminoAcidPos[j][l+1]):
                            return True
        
        return False
        
        
    def saveBestStructureUnits(self, outputFile, bestUnitList):
        """Ulozi nalezene nejvetsi vzajemne nekolidujici strukturni celky do vystupniho souboru."""

        # Zadny strukturni element nenalezen - do souboru se pouze zapise "nic"
        if(len(bestUnitList) == 0):
            outputFile = file(outputFile, "w")
            outputFile.write("")
            return
            
        # Sestaveni zahlavi - jmena proteinu
        lines = ""
        for i in range(len(bestUnitList[0].aminoAcidPos)):
            lines += bestUnitList[0].aminoAcidPos[i][0] + " "
        lines += "\n"
        
        # Sestaveni hlavniho obsahu - pozice aminokyselin nejvetsich navzajem nekolidujicich strukturnich celku
        for i in range(len(bestUnitList)):
            for j in range(len(bestUnitList[i].aminoAcidPos[0]) - 1):
                lines += bestUnitList[i].aminoAcidSeq[j] + " "
                for k in range(len(bestUnitList[i].aminoAcidPos)):
                    lines += str(bestUnitList[i].aminoAcidPos[k][j+1]) + " " 
                lines += "\n"

        # Ulozeni nejvetsich navzajem nekolidujicich strukturnich celku do vystupniho souboru
        outputFile = file(outputFile, "w")
        outputFile.write(lines)
        
    def showBestStructureUnits(self, bestUnitList, proteinId, proteinPath):
        """Zobrazi ziskane strukturni celky na vybranem proteinu."""
         
        # Nacteni souboru s proteinem
        cmd.load(proteinPath, "protein")
        cmd.set("sphere_scale", 0.3);
        
        # Vypocet oznacenych aminokyselin na proteinu
        aminoAcidList = []
        chain = "A"
        for i in range(len(bestUnitList)):
            aminoAcidUnit = []
            for j in range(len(bestUnitList[i].aminoAcidPos)):
                # Nalezeni proteinu mezi seznamy pozic
                if(bestUnitList[i].aminoAcidPos[j][0].find(proteinId) != -1):
                    pos = bestUnitList[i].aminoAcidPos[j][0].find("_")
                    chain = "A"
                    if(pos != -1):
                        chain = bestUnitList[i].aminoAcidPos[j][0][(pos+1):(pos+2)]
                        
                    for k in range(len(bestUnitList[i].aminoAcidPos[j]) - 1):
                        #print "///" + chain + "/" + str(bestUnitList[i].aminoAcidPos[j][k+1])
                        aminoAcidUnit.append("///" + chain + "/" + str(bestUnitList[i].aminoAcidPos[j][k+1]))
            aminoAcidList.append(aminoAcidUnit)

        # Zobrazeni strukturnich elementu
        color = ["red", "green", "yellow", "magenta", "cyan", "blue"]
        cmd.color("gray",  "all")
        cmd.color("white",  ("///" + chain))
        cmd.show("spheres",  "all") 
        
        for i in range(len(aminoAcidList)):
            for aminoAcid in aminoAcidList[i]:
                actIndex = i % len(color)
                actColor = color[i % len(color)]
                #cmd.show("spheres", aminoAcid)
                #cmd.show("spheres", aminoAcid)
                #cmd.show("ribbon", aminoAcid + "/CA")
                #cmd.color(actColor, aminoAcid)       
                cmd.color("red", aminoAcid)       
        # Zoomovani na finalni zobrazeni
        cmd.orient()

        
        
