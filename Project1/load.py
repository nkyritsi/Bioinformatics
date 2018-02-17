
###Takes fasta file name as input and returns a single string of all the DNA in
###the file.
def loadSeq(filename):
    infile = open(filename, "r")
    #print(infile.read())

    infile.readline() #get rid of title

    allLines = infile.read()
    allLines = allLines.replace('\n', '')
    return allLines


