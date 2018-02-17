import sys

#Function to load in fasta files
def loadSeq(filename):
    infile = open(filename, "r")

    title = infile.readline()
    allLines = infile.read()
    allLines = allLines.replace('\n', '')
    return (title, allLines)

def makeMatrices(fileName):
    matrixFile = open(fileName, "r")
    matrixFile.readline()
    matrixFile.readline() #ignore first 2 lines
    proteins_dna = matrixFile.readline()

    proteins_dna.split(" ")
    p = list(proteins_dna)
    while " " in p:
        p.remove(" ")

    Matrix = []
    Matrix.append(p)

    for line in matrixFile:
        l = list(line.split(' '))
        l.remove(l[0])
        while '' in l:
            l.remove('')
        while '\n' in l:
            l.remove('\n')
        for i in range(len(l)):
            l[i] = int(l[i])
        Matrix.append(l)

    return Matrix

#returns tuple of (seq1, seq1 positions, seq2, seq2 positions, alignment score)
def globalAlignment(seq1, seq2, gappenalty, scoreMatrix):

    n = len(seq1) + 1 #num rows/long
    m = len(seq2) + 1 #num cols/wide
    seqMatrix = [0] * n
    for i in range(n):
            seqMatrix[i] = [0] * m

    for i in range(1,n): #initialize matrix
        for j in range(1,m):
            if i-1 == 0:
                seqMatrix[0][j] = -j
            if j-1 == 0:
                seqMatrix[i][0] = -i

    seq1Bases = [0]+list(seq1) #to keep track of positions of bases
    seq2Bases = [0]+list(seq2)

    pathDict = {} #dictionary of coordinate keys with values as previous coordinates (to make path)
    prevCoord = (0, 0)
    for i in range(1,n):
        for j in range(1,m):

            base1 = seq1Bases[i]
            base2 = seq2Bases[j]
            base1index = scoreMatrix[0].index(base1)
            base2index = scoreMatrix[0].index(base2)
            matchscore = scoreMatrix[base1index+1][base2index] #pull score from matrix

            #three positions to look at
            upLeft = seqMatrix[i-1][j-1] + matchscore
            up = seqMatrix[i-1][j] + gappenalty
            left = seqMatrix[i][j-1] + gappenalty

            if max(upLeft, up, left) == upLeft:
                seqMatrix[i][j] = upLeft
                prevCoord = (i-1, j-1)

            if max(upLeft, up, left) == up:
                seqMatrix[i][j] = up
                prevCoord = (i-1, j)

            if max(upLeft, up, left) == left:
                seqMatrix[i][j] = left
                prevCoord = (i, j-1)

            pathDict[(i, j)] = prevCoord

    np = n-1
    mp = m-1
    finalCoord = (np, mp)
    finalScore = seqMatrix[np][mp]
    s1 = [] #lists to add bases to (for returning alignment)
    s2 = []

    s1EndPos = finalCoord[0]
    s2EndPos = finalCoord[1]

    while finalCoord in pathDict:
        currPos = pathDict.get(finalCoord)
        if currPos == (finalCoord[0]-1, finalCoord[1]-1):#if prev is upper left
            s1.append(seq1Bases[finalCoord[0]])
            s2.append(seq2Bases[finalCoord[1]])
            finalCoord = currPos
        elif currPos == (finalCoord[0]-1, finalCoord[1]):
            s1.append(seq1Bases[finalCoord[0]-1])
            s2.append("-")
            finalCoord = currPos
        elif currPos == (finalCoord[0], finalCoord[1]-1):
            s1.append("-")
            s2.append(seq2Bases[finalCoord[1]-1])
            finalCoord = currPos

    s1StartPos = finalCoord[0]
    s2StartPos = finalCoord[1]

    s1Positions = (s1StartPos, s1EndPos)
    s2Positions = (s2StartPos, s2EndPos)

    s1.reverse()
    s2.reverse()

    return (s1, s1Positions, s2, s2Positions, finalScore)

def semiGlobalAlignment(seq1, seq2, gappenalty, scoreMatrix):
    n = len(seq1) + 1 #num rows/long
    m = len(seq2) + 1 #num cols/wide
    seqMatrix = [0] * n
    for i in range(n):
            seqMatrix[i] = [0] * m

    seq1Bases = [0]+list(seq1) #to keep track of positions of bases
    seq2Bases = [0]+list(seq2)

    pathDict = {} #dictionary of coordinate keys with values as previous coordinates (to make path)
    prevCoord = (0, 0)

    #add scores to seqMatrix
    for i in range(1,n):
        for j in range(1,m):

            base1 = seq1Bases[i]
            base2 = seq2Bases[j]
            base1index = scoreMatrix[0].index(base1)
            base2index = scoreMatrix[0].index(base2)
            matchscore = scoreMatrix[base1index+1][base2index] #pull score from matrix

            #three positions to look at
            upLeft = seqMatrix[i-1][j-1] + matchscore
            up = seqMatrix[i-1][j] + gappenalty
            left = seqMatrix[i][j-1] + gappenalty

            #if in last column or last row, no gap penalty
            if i == n-1 or j == m:
                up = seqMatrix[i-1][j]
                left = seqMatrix[i][j-1]

            if max(upLeft, up, left) == upLeft:
                seqMatrix[i][j] = upLeft
                prevCoord = (i-1, j-1)

            if max(upLeft, up, left) == up:
                seqMatrix[i][j] = up
                prevCoord = (i-1, j)

            if max(upLeft, up, left) == left:
                seqMatrix[i][j] = left
                prevCoord = (i, j-1)

            pathDict[(i, j)] = prevCoord

    np = n-1
    mp = m-1
    finalCoord = (np, mp)
    finalScore = seqMatrix[np][mp]
    s1 = [] #lists to add bases to (for returning alignment)
    s2 = []

    s1EndPos = finalCoord[0]
    s2EndPos = finalCoord[1]

    while finalCoord in pathDict:
        currPos = pathDict.get(finalCoord)
        if currPos == (finalCoord[0]-1, finalCoord[1]-1):#if prev is upper left
            s1.append(seq1Bases[finalCoord[0]])
            s2.append(seq2Bases[finalCoord[1]])
            finalCoord = currPos
        elif currPos == (finalCoord[0]-1, finalCoord[1]):
            s1.append(seq1Bases[finalCoord[0]-1])
            s2.append("-")
            finalCoord = currPos
        elif currPos == (finalCoord[0], finalCoord[1]-1):
            s1.append("-")
            s2.append(seq2Bases[finalCoord[1]-1])
            finalCoord = currPos

    s1StartPos = finalCoord[0]
    s2StartPos = finalCoord[1]

    s1Positions = (s1StartPos, s1EndPos)
    s2Positions = (s2StartPos, s2EndPos)

    s1.reverse()
    s2.reverse()

    return (s1, s1Positions, s2, s2Positions, finalScore)

#returns tuple of first sequence, first eq positions, second sequence, second seq positions, and the alignment score
def localAlignment(seq1, seq2, gappenalty, scoreMatrix):
    n = len(seq1) + 1 #num rows/long
    m = len(seq2) + 1 #num cols/wide
    seqMatrix = [0] * n
    for i in range(n):
            seqMatrix[i] = [0] * m

    seq1Bases = [0]+list(seq1) #to keep track of positions of bases
    seq2Bases = [0]+list(seq2)

    pathDict = {} #dictionary of coordinate keys with values as previous coordinates (to make path)
    prevCoord = (0, 0)

    maxScore = 0
    maxScorePositions = (0,0)

    #add scores to seqMatrix
    for i in range(1,n):
        for j in range(1,m):

            base1 = seq1Bases[i]
            base2 = seq2Bases[j]
            base1index = scoreMatrix[0].index(base1)
            base2index = scoreMatrix[0].index(base2)
            matchscore = scoreMatrix[base1index+1][base2index] #pull score from matrix

            #three positions to look at
            upLeft = seqMatrix[i-1][j-1] + matchscore
            up = seqMatrix[i-1][j] + gappenalty
            left = seqMatrix[i][j-1] + gappenalty

            #if no positive values, score for position is 0, update max score to find where to start backtracking
            if max(upLeft, up, left) == upLeft:
                if upLeft >= 0:
                    seqMatrix[i][j] = upLeft
                else:
                    seqMatrix[i][j] = 0
                prevCoord = (i-1, j-1)

                if upLeft > maxScore:
                    maxScore = upLeft
                    maxScorePositions = (i, j)

            if max(upLeft, up, left) == up:
                if up >= 0:
                    seqMatrix[i][j] = up
                else:
                    seqMatrix[i][j] = 0
                prevCoord = (i-1, j)

                if up > maxScore:
                    maxScore = up
                    maxScorePositions = (i, j)

            if max(upLeft, up, left) == left:
                if left >= 0:
                    seqMatrix[i][j] = left
                else:
                    seqMatrix[i][j] = 0
                prevCoord = (i, j-1)

                if left > maxScore:
                    maxScore = left
                    maxScorePositions = (i, j)

            pathDict[(i, j)] = prevCoord

    finalCoord = maxScorePositions
    s1 = [] #lists to add bases to (for returning alignment)
    s2 = []

    s1EndPos = finalCoord[0]
    s2EndPos = finalCoord[1]

    while finalCoord in pathDict and seqMatrix[finalCoord[0]][finalCoord[1]] != 0: #stop backtracking when score for square is 0
        currPos = pathDict.get(finalCoord)
        if currPos == (finalCoord[0]-1, finalCoord[1]-1):#if prev is upper left
            s1.append(seq1Bases[finalCoord[0]])
            s2.append(seq2Bases[finalCoord[1]])
            finalCoord = currPos
        elif currPos == (finalCoord[0]-1, finalCoord[1]):
            s1.append(seq1Bases[finalCoord[0]-1])
            s2.append("-")
            finalCoord = currPos
        elif currPos == (finalCoord[0], finalCoord[1]-1):
            s1.append("-")
            s2.append(seq2Bases[finalCoord[1]-1])
            finalCoord = currPos

    s1StartPos = finalCoord[0]
    s2StartPos = finalCoord[1]

    s1Positions = (s1StartPos, s1EndPos)
    s2Positions = (s2StartPos, s2EndPos)

    s1.reverse()
    s2.reverse()

    return (s1, s1Positions, s2, s2Positions, maxScore)



def main():

    firstSeqIdx = sys.argv.index("-i") + 1
    firstSeqFile = str(sys.argv[firstSeqIdx])
    secSeqIdx = sys.argv.index("-j") + 1
    secSeqFile = str(sys.argv[secSeqIdx])

    seq1 = loadSeq(firstSeqFile) #first part of tuple is title name
    seq2 = loadSeq(secSeqFile)  #first part of tuple is title name

    seqType = sys.argv.index("-p") + 1 #T or F
    alignType = sys.argv.index("-atype") + 1 

    aType = sys.argv[alignType]
    sType = sys.argv[seqType]

    outFile = sys.argv.index("-o") + 1

    scorematrix = 0
    penalty = 0
    if sType == "T":
        scorematrix = makeMatrices("BLOSUM45")
        penalty = -5
    elif sType == "F":
        scorematrix = makeMatrices("dnaMatrix")
        penalty = -10

    alignment = ''
    if aType == "G":
        alignment = globalAlignment(seq1[1], seq2[1], penalty, scorematrix) #(seq1, seq1pos, seq2, seq2pos, score)
    elif aType == "S":
        alignment = semiGlobalAlignment(seq1[1], seq2[1], penalty, scorematrix) #(seq1, seq1pos, seq2, seq2pos, score)
    elif aType == "L":
        alignment = localAlignment(seq1[1], seq2[1], penalty, scorematrix) #(seq1, seq1pos, seq2, seq2pos, score)

    sameBases = 0
    for i in range(len(alignment[0])):
        if alignment[0][i] == alignment[2][i] and alignment[2][i] != "-":
            sameBases += 1
    identity = str(sameBases) + "/" + str(len(alignment[0]))

    seq1Title = seq1[0]
    name1 = seq1Title.split(" ")
    seq2Title = seq2[0]
    name2 = seq2Title.split(" ")

    a1 = ''.join(alignment[0])
    a2 = ''.join(alignment[2])
    alignmentFile = open(sys.argv[outFile], "w")

    if len(a1) > 60:
        rangeFile = len(a1)//60
        alignmentFile.write(str(name1[0][:-1]) + ": " + str(alignment[1][0]) + " " + a1[:61] + " " + str(alignment[1][0] + 60) + "\n")
        alignmentFile.write(str(name1[0][:-1]) + ": " + str(alignment[3][0]) + " " + a2[:61] + " " + str(alignment[3][0] + 60) + "\n" + "\n")
        j=60
        for i in range(rangeFile):
            if i != rangeFile-1:
                alignmentFile.write(name1[0][:-1] + ": " + str(alignment[1][0] + (j+1)) + " " + a1[j:j+61] + " " + str(alignment[1][0] + (j+1)+60) + "\n")
                alignmentFile.write(name2[0][:-1] + ": " + str(alignment[3][0] + (j+1)) + " " + a2[j:j+61] + " " + str(alignment[3][0] + (j+1)+60) + "\n" + "\n")
                j+=61
            else:
                alignmentFile.write(name1[0][:-1] + ": " + str(alignment[1][0] + (j+1)) + " " + a1[j:j+61] + " " + str(alignment[1][1]) + "\n")
                alignmentFile.write(name2[0][:-1] + ": " + str(alignment[3][0] + (j+1)) + " " + a2[j:j+61] + " " + str(alignment[3][1]) + "\n" + "\n")

        alignmentFile.write("Score: " + str(alignment[4]) + "\n")
        alignmentFile.write("Identity: " + str(identity) + "\n")
    else:
        alignmentFile.write(name1[0][:-1] + " " + str(alignment[1][0]) + " " + a1 + " " + str(alignment[1][1]) + "\n")
        alignmentFile.write(name2[0][:-1] + " " + str(alignment[3][0]) + " " + a2 + " " + str(alignment[3][1]) + "\n")
        alignmentFile.write("Score: " + str(alignment[4]) + "\n")
        alignmentFile.write("Identity: " + str(identity) + "\n")

main()














