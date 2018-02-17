from load import *
from dna import *
import random


##Starts at 0 position of DNA and searches in units of 3 looking
##for start codons. When it finds a start codon, take the slice
##of DNA beginning with "ATG" until a stop codon. Store the sequence
##in a list and continue searching for next "ATG". Returns all
##ORFs found.
def oneFrame(DNA):
    listORF = []

    for i in range(len(DNA)):
        if DNA[i:i+3] == 'ATG':
            newORF = ''

            j = i
            while (DNA[j:j+3] != 'TAG'and DNA[j:j+3] != 'TGA' and DNA[j:j+3] != 'TAA') and j < len(DNA):
                #print(DNA[j])
                newORF += DNA[j:j+3]
                j += 3

            listORF.append(newORF)
            newORF = ''

    return listORF


def oneFrameV2(DNA):
    listORF = []

    for i in range(len(DNA)):
        if DNA[i:i+3] == 'ATG':
            newORF = ''

            j = i
            while (DNA[j:j+3] != 'TAG'and DNA[j:j+3] != 'TGA' and DNA[j:j+3] != 'TAA') and j < len(DNA):
                #print(DNA[j])
                newORF += DNA[j:j+3]
                j += 3

            listORF.append(newORF)
            if len(listORF) >1 and newORF in listORF[-2]:
                listORF.remove(newORF)
            newORF = ''

    return listORF

#Calls oneFrameV2 to find the longest open reading frame of all 3 reading frames
def longestORF(DNA):
    allORFS = oneFrameV2(DNA)
    longORF = max(allORFS, key=len)
    return longORF

#Finds the longest ORF on DNA and its reverse complement
def longestORFBothStrands(DNA):
    longestPositiveStrand = longestORF(DNA)

    reverseStrand = reverseComplement(DNA)
    longestReverseStrand = longestORF(reverseStrand)

    if len(longestPositiveStrand) >= len(longestReverseStrand):
        return longestPositiveStrand
    else:
        return longestReverseStrand

def collapse(DNAlist):
    DNAstring = ''

    for i in DNAlist:
        DNAstring += i

    return DNAstring

#Makes garbage sequences, finds the longest garbage ORF, and returns its length
def longestORFNoncoding(DNA, numReps):
    DNAlist = list(DNA)
    ORFsizes = []

    for i in range(numReps):
        random.shuffle(DNAlist)
        DNAs = collapse(DNAlist)

        longORF = longestORFBothStrands(DNAs)
        ORFsizes.append(len(longORF))

    return max(ORFsizes)

#Identify all ORFs in real (un-shuffled) DNA and return them as a list.
##If there are none, return an empty list
def findORFs(DNA):
    allORFS = []

    ORFSframe1 = oneFrameV2(DNA)
    #ORFSframe2 = oneFrameV2(DNA[1:])
    #ORFSframe3 = oneFrameV2(DNA[2:])

    allORFS = allORFS + ORFSframe1 #+ ORFSframe2 + ORFSframe3
    return allORFS

#Searches both forward and reverse strands for ORFs and returns a list with
#all the ORFs found
def findORFsBothStrands(DNA):
    allORFS = []

    ORFsForward = oneFrameV2(DNA)

    reverseStrand = reverseComplement(DNA)
    ORFsBackward = oneFrameV2(reverseStrand)

    allORFS += (ORFsBackward + ORFsForward)
    return allORFS

#returns the beginning and end coordinates of an ORF in DNA
def getCoordinates(orf, DNA):
    ORFlength = len(orf)
    beginningCoord = DNA.find(orf)

    if beginningCoord == -1:
        reverseORF = reverseComplement(orf)
        beginningCoord = DNA.find(reverseORF)

    endCoord = beginningCoord + ORFlength
    return [beginningCoord, endCoord]

#Identifies ORFs longer than minLen, and returns list of
# [beginningCoord, endCoord, proteinSeq]
def geneFinder(DNA, minLen):
    longerORFs = []

    allORFS = findORFsBothStrands(DNA)
    for i in allORFS:
        if len(i) > minLen:
            longerORFs.append(i)

    finalOutputList = []
    for orf in longerORFs:
        coords = getCoordinates(orf, DNA)
        proteinSeq = codingStrandToAA(orf)
        orfInfo = [coords[0], coords[1], proteinSeq]
        finalOutputList.append(orfInfo)

    finalOutputList.sort()
    return finalOutputList

def printGenes(geneList):
    for i in geneList:
        print(i)

