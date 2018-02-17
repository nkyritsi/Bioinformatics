##contains helper functions for manipulating DNA sequences

from aminoAcids import *

##Takes DNA string as input and returns the sequence of its complementary
##DNA strand, also in 5' to 3' order
def reverseComplement(DNA):
    reverseString = ''
    for i in range(1,len(DNA)+1):
        lastBase = DNA[-i]
        if lastBase == 'A':
            reverseString += 'T'
        elif lastBase == 'T':
            reverseString += 'A'
        elif lastBase == 'C':
            reverseString += 'G'
        else:
            reverseString += 'C'
    return reverseString

##Takes a DNA sequence and returns the corresponding amino acids as a string.
##If the length of DNA sequence is not divisible by 3, print an error and return NONE
def codingStrandToAA(DNA):

    stringAA = ''

    DNAlength = len(DNA)
    if (DNAlength%3) != 0:
        print("DNA sequence is not divisible by 3!")
        return None
    else:
        rangeLoop = DNAlength//3 #Number of codons to iterate through
        j = 0
        for i in range(rangeLoop):
            codon = DNA[j:j+3]
            j +=3

            #See which index codon is in
            for k in range(len(codons)):
                if codon in codons[k]:
                    stringAA += aa[k]
    return stringAA