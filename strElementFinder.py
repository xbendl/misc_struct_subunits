#!/usr/bin/python
# -*- coding: utf-8 -*-

# PBI - Projekt, varianta B
#
# strElementFinder.py - Modul pro identifikaci strukturnich elementu - dvou
#                       paru aminokyselin sousedicich na kostre, pricemz 
#                       vzdalenost mezi proteiny je mensi nez 5 angstromu.
#
# Autor: Jaroslav Bendl, xbendl00
#

from aminoAcidRecord import AminoAcidRecord
from strElement import StrElement

class StrElementFinder:
    """Ziska strukturni elementy shodne pro zadanou sadu proteinu."""
    
    def __init__(self, angstromDist):
        self.elementList = []             # list pro ulozeni nalezenych strukturnich elementu
        self.angstromDist = angstromDist  # vzdalenost v angstromech, kdy se dva c-alpha uhliky na kostre berou jako sousedni

    def __getAminoList(self, proteinStructure):
        """Ziska seznam aminokyselin s pozici C-alpha uhliku pro danou proteinovou strukturu."""
        actAminoList = []
        infile = open(proteinStructure.proteinFile, "r")
        modelFound = False
        
        # Ziskani posloupnosti aminokyselin na retezci s informaci o pozici C-alpha uhliku
        while infile:
            line = infile.readline()
            if(len(line) == 0):            # Nic nenacteno => konec souboru
                break 

            # Zabraneni zpracovani vice nez jednoho modelu (bere se pouze ten v poradi prvni)
            if((line[0:5] == "MODEL") and (modelFound == True)):
                # Kontrola zpracovani pozadovaneho retezce na strukture proteinu
                break
            
            if((line[0:5] == "MODEL")):
                modelFound = True
                
            # Zpracovani radku s oznacenim ATOM (prochazeni po pateri)
            if(line[0:4] == "ATOM"):
                # Kontrola zpracovani pozadovaneho retezce na strukture proteinu
                if((len(proteinStructure.proteinChain) == 1) and (proteinStructure.proteinChain != line[21:22])):
                    continue
                   
                # Pro kazdou aminokyselinu se pridava zaznam jen o jednom atomu - C-alpha uhliku
                if(line[12:16].strip() == "CA"):
                    actAminoList.append(AminoAcidRecord(int(line[23:26]), line[17:20], float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()), self.angstromDist))

        return actAminoList


    def getStructureElements(self, proteinList):
        """Ziska strukturni elementy spolecne pro vstupni sadu proteinu."""
        
        self.elementList = []                  # list pro ulozeni nalezenych strukturnich elementu

        # Kontrola minimalniho poctu proteinovych struktur
        if(len(proteinList) < 2):
            print "CHYBA! Na vstupu musi byt alespon dve korektni proteinove struktury."
            exit()
                
        # Ziskani seznamu posloupnosti aminokyselin s pozici C-alpha uhliku pro kazdou vysetrovanou proteinovou strukturu
        aminoList = []
        for i in range(len(proteinList)):
            aminoList.append(self.__getAminoList(proteinList[i]))
            
            # Nastaveni vzdalenosti k sousednimu C-alpha uhliku (pro urychleni vypoctu)
            for j in range(len(aminoList[i]) - 1):
                aminoList[i][j].setNeighbourDistance(aminoList[i][j+1])
        
        
        # Vyhledavani prvniho paru z posloupnosti prvni proteinove struktury
        for i in range(len(aminoList[0]) - 3):
            i1 = i
            i2 = i + 1
            
            # Vyhledavani druheho paru z posloupnosti prvni proteinove struktury
            for j in range(len(aminoList[0]) - 3 - i):
                j1 = i + j + 2
                j2 = i + j + 3
                
                # Kontrola vzdalenosti mezi c-alpha u i1-i2, j1-j2, i1-j1, i1-j2, i2-j1, i2-j2
                if((aminoList[0][i1].isNearNeighbour() == True) and (aminoList[0][j1].isNearNeighbour() == True) and \
                (aminoList[0][i1].isNear(aminoList[0][j1])) and (aminoList[0][i1].isNear(aminoList[0][j2])) and \
                (aminoList[0][i2].isNear(aminoList[0][j1])) and (aminoList[0][i2].isNear(aminoList[0][j2]))):
                    
                    # Vytvoreni noveho strukturniho elementu a pridani pozic vyskytu pro prvni proteinovou strukturu
                    element = StrElement(aminoList[0][i1].aminoAcid + aminoList[0][i2].aminoAcid, aminoList[0][j1].aminoAcid + aminoList[0][j2].aminoAcid, [])
                    proteinName = proteinList[0].proteinId
                    if(len(proteinList[0].proteinChain) == 1):
                        proteinName = proteinName + "_" + proteinList[0].proteinChain
                    element.addProteinPosition(proteinName, aminoList[0][i1].resiNumber, aminoList[0][i2].resiNumber, aminoList[0][j1].resiNumber, aminoList[0][j2].resiNumber)
                    
                    #print "Kandidat: " + aminoList[0][i1].aminoAcid + aminoList[0][i2].aminoAcid + " " + aminoList[0][j1].aminoAcid + aminoList[0][j2].aminoAcid + "|" + str(aminoList[0][i1].resiNumber) + ":" + str(aminoList[0][i2].resiNumber) + ":" + str(aminoList[0][j1].resiNumber) + ":" + str(aminoList[0][j2].resiNumber)
                    
                    # Vyhledavani v dalsich proteinovych strukturach
                    found = False
                    for k in range(len(proteinList) - 1):
                        
                        found = False
                        # Vyhledavani prvniho paru druhe proteinove struktury
                        for l in range(len(aminoList[k+1]) - 3):                            
                            l1 = l
                            l2 = l + 1
                            
                            # Kontrola shodnosti prvniho paru v druhe proteinove strukture
                            if((aminoList[k+1][l1].aminoAcid == aminoList[0][i1].aminoAcid) and \
                                (aminoList[k+1][l2].aminoAcid == aminoList[0][i2].aminoAcid)):
                                
                                for m in range(len(aminoList[k+1]) - 3 - l):
                                    m1 = l + m + 2
                                    m2 = l + m + 3
                                    
                                    # Kontrola shodnosti druheho paru v druhe proteinove strukture
                                    if((aminoList[k+1][m1].aminoAcid == aminoList[0][j1].aminoAcid) and \
                                        (aminoList[k+1][m2].aminoAcid == aminoList[0][j2].aminoAcid)):
                                                                              
                                        # Kontrola vzdalenosti mezi c-alpha u l1-l2, m1-m2, l1-m1, l1-m2, l2-m1, l2-m2
                                        if((aminoList[k+1][l1].isNearNeighbour() == True)  and \
                                            (aminoList[k+1][m1].isNearNeighbour() == True) and \
                                            (aminoList[k+1][l1].isNear(aminoList[k+1][m1]))  and \
                                            (aminoList[k+1][l1].isNear(aminoList[k+1][m2]))  and \
                                            (aminoList[k+1][l2].isNear(aminoList[k+1][m1]))  and \
                                            (aminoList[k+1][l2].isNear(aminoList[k+1][m2]))):
                                                
                                            # Nastaveni priznaku existence strukturniho elementu 
                                            found = True
                                            
                                            # Pridani pozic vyskutu elementu pro tuto proteinovou strukturu
                                            proteinName = proteinList[k+1].proteinId
                                            if(len(proteinList[k+1].proteinChain) == 1):
                                                proteinName = proteinName + "_" + proteinList[k+1].proteinChain
                                            element.addProteinPosition(proteinName, aminoList[k+1][l1].resiNumber, aminoList[k+1][l2].resiNumber, aminoList[k+1][m1].resiNumber, aminoList[k+1][m2].resiNumber)
                                            
                        # Element nebyl nalezen v posledni zkomane proteinove strukture, neni tedy platny
                        if(found == False):
                            break
                    
                    # Pridani elementu do seznamu elementu (element byl nalezen ve vsech proteinovych strukturach)
                    if(found == True):
                        self.elementList.append(element)

        # Uprava seznamu elementu - rozepsani vicenasobnych vyskytu druheho paru aminokyselin v ostatnich strukturach
        i = 0
        while(i < len(self.elementList)):
            lastProteinId = "X" # pomocna promenna pro detekci vice vyskytu druheho paru aminokyselin v jedne strukture
            multipleOccurence = True # pomocna promenna pro signalizace problemu s vicenasobnym vyskytem
            
            # Cyklus fixujici pozici - zarucujici kompletni rozbaleni vicenasobneho vyskytu
            while(multipleOccurence == True):
                multipleOccurence = False
                
                # Pruchod pres vsechny vyskyty 
                for j in range(len(self.elementList[i].proteinPosition)):
                    # Detekce opakovaneho vyskytu druheho paru aminokyselin v jedne strukture
                    if(lastProteinId == self.elementList[i].proteinPosition[j][0]):
                        # Nalezeni koncoveho indexu sekvence vicenasobneho vyskytu
                        k = 0
                        
                        while k < (len(self.elementList[i].proteinPosition) - j):
                            if(self.elementList[i].proteinPosition[j+k][0] != lastProteinId):
                                break;
                            k = k+1
                        
                        # Vytvoreni dvou nahrazujicich sekvenci vyskytu, jejich umisteni a vyjmuti te puvodni
                        childElement1 = StrElement(self.elementList[i].aminoAcidPair1, self.elementList[i].aminoAcidPair2, (self.elementList[i].proteinPosition[0:j] + self.elementList[i].proteinPosition[j+k:]))
                        
                        childElement2 = StrElement(self.elementList[i].aminoAcidPair1, self.elementList[i].aminoAcidPair2, (self.elementList[i].proteinPosition[0:(j-1)] + self.elementList[i].proteinPosition[j:]))
                        
                        self.elementList.pop(i)
                        self.elementList.insert(i, childElement2)
                        self.elementList.insert(i, childElement1)
                        
                        multipleOccurence = True  # nastaveni priznaku - je nutne na danem indexu jeste zustat
                        break;
                        
                    lastProteinId = self.elementList[i].proteinPosition[j][0]
                    
            i = i + 1
            
        return self.elementList


    def saveElements(self, outputFile):
        """Ulozi nalezene strukturni elementy do vystupniho souboru."""
    
        # Ziskani strukturnich elementu
        lines = ""
        for i in range(len(self.elementList)):
            lines += self.elementList[i].getStructureElement() + "\n"
                
        # Ulozeni strukturnich elementu do vystupniho souboru
        outputFile = file(outputFile, "w")
        outputFile.write(lines);


