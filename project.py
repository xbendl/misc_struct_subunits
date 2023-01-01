#!/usr/bin/python
# -*- coding: utf-8 -*-

# PBI - Projekt, varianta B
#
# project.py - Identifikace spolecnych strukturnich elementu vysetrovanych 
#              proteinu, jejich spojeni do vetsich celku a jejich zobrazeni na
#              libovolnem proteinu.
#
# Autor: Jaroslav Bendl, xbendl00
#

import sys

from parameter import Parameter
from record import Record
from strElementFinder import StrElementFinder
from strElementJoiner import StrElementJoiner


def main():
    """Ridi cinnost programu."""
    
    # Nacte seznam vysetrovanych proteinu a zajisti jejich dostupnost ve slozce proteins
    param = Parameter("configuration.txt")
    proteinList = param.getProteinList()

    # Ziskani strukturnich elementu a jejich ulozeni do souboru
    strElemFinder = StrElementFinder(param.angstromDist)
    strElemList   = strElemFinder.getStructureElements(proteinList)
    print "> Identifikace spolecnych strukurnich elementu probehla uspesne."
    strElemFinder.saveElements(param.outputFile1)
    print "> Seznam spolecnych strukturnich elementu byl ulozen do souboru: " + param.outputFile1

    # Ziskani strukturnich celku a jejich ulozeni do souboru
    strElemJoiner = StrElementJoiner()
    strElemJoiner.getStructureUnits(strElemList)
    bestUnitList = strElemJoiner.getBestStructureUnits(param.unitCount)
    print "> Identifikace strukurnich celku probehla uspesne."
    strElemJoiner.saveBestStructureUnits(param.outputFile2, bestUnitList)
    proteinPath = param.getShowProteinPath()
    strElemJoiner.showBestStructureUnits(bestUnitList, param.showProteinId, proteinPath)
    print "> Seznam nejvesich strukturnich celku byl ulozen do souboru: " + param.outputFile2

main()   # GO!
