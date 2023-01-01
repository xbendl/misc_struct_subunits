#!/usr/bin/python
# -*- coding: utf-8 -*-

# PBI - Projekt, varianta B
#
# parameter.py - Modul pro cteni vstupniho souboru s identifikatory
#                vysetrovanych proteinovych struktur a jejich pripadne 
#                stazeni z internetu.
#
# Autor: Jaroslav Bendl, xbendl00
#

from httplib import HTTPConnection
from ftplib import FTP
from record import Record

import sys
import gzip
import os.path

class Parameter:
    """Trida zajistujici dostupnost pozadovanych struktur z databaze PDB."""
    
    def __init__(self, inputFile):
        self.outputFile1   = ""  # cesta k vystupnimu souboru pro ulozeni strukturnich elementu
        self.outputFile2   = ""  # cesta k vystupnimu souboru pro ulozeni strukturnich celku
        self.showProteinId = ""  # identifikator proteinu, na kterem budou vyznaceny nalezene strukturni celky
        self.unitCount     = 0   # pocet uvazovanych nejvetsich navzajem nekolidujicich strukturnich struktur
        self.angstromDist  = 5   # vzdalenost v angstromech, kdy se dva c-alpha uhliky na kostre berou jako sousedni
        self.proteinList   = []  # seznam s identifikatory zkoumanych proteinu
       
        self.httpAddress = "dx.doi.org"
        self.fileHttpPrefix = "/10.2210/pdb"
        self.fileHttpSuffix = "/pdb"
        self.proteinDirectory = "./proteins/"
        self.proteinFilePrefix = "pdb"
        self.proteinFileSuffix = ".ent"

        self.__readInputFile(inputFile)


    def __download(self, proteinId):
        """Stahne protein z databaze PDB."""
        
        # Stahnuti webove stranky s korektni adresou zaznamu PDB na FTP serveru
        conn = HTTPConnection(self.httpAddress)
        conn.request("GET", (self.fileHttpPrefix + proteinId + self.fileHttpSuffix))
        r = conn.getresponse()
        data = r.read()
        conn.close()
                
        # Osetreni situace, ze zaznam z danym identifikatorem neni v databazi PDB
        if(data.find("Error - DOI Not Found") != -1):
            print "VAROVANI! Zaznam s identifikatorem " + proteinId + " neni v databazi PDB."
            return

        # Vytahnuti korektni adresy z nactene webove stranky k pozadovanemu PDB zaznamu
        startFtpPos = data.find("=") + 2        
        endFtpPos = data.find("\"", startFtpPos)
        ftpAddressFull = data[startFtpPos:endFtpPos]        

        # Stahnuti archivu se zaznamem z databaze PDB na FTP serveru
        ftpAddressFull = ftpAddressFull[len("ftp://"):]
        ftpAddress = ftpAddressFull[0 : ftpAddressFull.find("/")]
        fileFtp = ftpAddressFull[len(ftpAddress) : ]
        
        fileFtpGz = proteinId + ".gz"
        ftpConn = ftp = FTP(ftpAddress)
        ftpConn.login()
        ftpConn.retrbinary("RETR " + fileFtp, open(self.proteinDirectory + fileFtpGz, "wb").write)
        ftpConn.quit()
        print "> Stazeni souboru se zaznamem " + proteinId + " probehlo uspesne."
        
        # Rozbaleni archivu se zaznamem a ulozeni nekomprimovane verze
        fileGz = gzip.open(self.proteinDirectory + fileFtpGz, "rb")
        data = fileGz.read()
        fileGz.close()
        
        outputFileName = self.proteinDirectory + self.proteinFilePrefix + proteinId + self.proteinFileSuffix
        outputFile = file(outputFileName, "w")
        outputFile.write(data);
        print "> Rozbaleni souboru se zaznamem " + proteinId + " probehlo uspesne."
        
        # Smazani archivu se zaznamem (dekomprimace jiz probehla, neni duvod jej dale uchovavat)
        os.remove(self.proteinDirectory + fileFtpGz)
        return


    def __readInputFile(self, inputFile):
        """Nacte z lokalniho umisteni pozadovane proteinove zaznamy, pripadne je stahne z internetu."""       
        infile = open(inputFile, "r")
        
        # Pruchod vstupnim souborem po radcich
        while infile:
            line = infile.readline()
            if(len(line) == 0):            # Nic nenacteno => konec souboru
                break
            if(len(line.strip()) == 0):    # Pouze prazdne znaky => preskocit na dalsi radek
                continue

            if(line[0:12] == "OUTPUT_FILE1"):    # Cteni radku s cestou k souboru pro ulozeni strukturnich jednotek
                self.outputFile1 = line[(line.find("=") + 1):].strip()
                continue
            if(line[0:12] == "OUTPUT_FILE2"):    # Cteni radku s cestou k souboru pro ulozeni strukturnich celku
                self.outputFile2 = line[(line.find("=") + 1):].strip()
                continue
            if(line[0:10]  == "UNIT_COUNT"):     # Cteni radku s poctem strukturnich celku uvazovanych pro vysledek
                self.unitCount = int(line[(line.find("=") + 1):].strip())
                continue
            if(line[0:15] == "SHOW_PROTEIN_ID"): # Cteni radku s identifikatorem proteinu pro vyznaceni struk. celku
                self.showProteinId = line[(line.find("=") + 1):].strip()
                continue
            if(line[0:13] == "DIST_ANGSTROM"):   # Cteni radku s identifikatorem proteinu pro vyznaceni struk. celku
                self.angstromDist = float(line[(line.find("=") + 1):].strip())
                continue
            if(line[0:12] == "PROTEIN_LIST"):    # Cteni radku uvozujiciho sekci s identifikatoru uvazovanych proteinu
                continue
           
            # Radky s identifikatory uvazovanych proteinu
            proteinId = line.strip()
            proteinChain = ''
            
            # Rozdeleni nacteneho radku na identifikator proteinu a oznaceni retezce
            sepPosition = proteinId.find("_")
            if(sepPosition != -1):
                proteinChain = proteinId[(sepPosition + 1) : (sepPosition + 2)].upper()
                proteinId = proteinId[0 : sepPosition]
                        
            # Kontrola spravneho formatu identifikatoru proteinu a oznaceni pripadneho retezce
            if((len(proteinId) != 4) or ((len(proteinChain) == 1) and (proteinChain.isalnum() == False))):
                print "VAROVANI! Nespravne zadany identifikator proteinu: " + line.strip()
                continue
            
            # Kontrola existence souboru s proteinem v prislusne slozce, pripadne jeho stazeni z internetu
            proteinPath = self.proteinDirectory + self.proteinFilePrefix + proteinId.lower() + self.proteinFileSuffix
            if(os.path.exists(proteinPath) == False):
                if(self.__download(proteinId.lower()) == False):
                    continue

            self.proteinList.append(Record(proteinId, proteinPath, proteinChain))
        
        # Kontrola korektnosti nastaveni parametru showProteinId z konfiguracniho souboru (musi se jednat o identifikator jednoho z vysetrovanych proteinu)
        found = False
        for i in range(len(self.proteinList)):
            proteinId = self.proteinList[i].proteinId
            if(proteinId.find("_") != -1):
                proteinId = proteinId[0:proteinId.find("_")]  # smazani pripadneho retezce v identifikatoru
            if(proteinId == self.showProteinId):
                found = True
                break
        if(found == False):
            print "CHYBA! Nespravne zadany parametr \"showProteinId\" v konfiguracnim souboru. Tento parametr musi mit hodnotu jednoho z vysetrovanych proteinu."
            sys.exit()
        
    def getProteinList(self):
        """Vrati seznam s vysetrovanymi proteiny."""
        return self.proteinList


    def getShowProteinPath(self):
        """Vrati cestu k PDB zaznamu proteinovu, na kterem se zobrazi nalezene nejvetsi strukturni celky."""
        showProteinPath = self.proteinDirectory + self.proteinFilePrefix + self.showProteinId.lower() + self.proteinFileSuffix
        return showProteinPath
        