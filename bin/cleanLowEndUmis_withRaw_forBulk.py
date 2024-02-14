import sys
import getopt
import os.path
import glob
import argparse

# Argument parsing:
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, metavar='<input UMI distributions>',
                        help='input UMI distributions file')
    parser.add_argument('-o', type=str, required=True, metavar='<output UMI counts file>',
                        help='output UMI counts file')
    parser.add_argument('-r', type=str, required=False, default=None, metavar='<output read counts file>',
                        help='output read counts file')
    parser.add_argument('-n', type=int, required=False, default=2, metavar='<minimum UMI count for merging>',
                        help='minimum required UMI count for merging singletons')
    parser.add_argument('-u', type=int, required=False, default=1, metavar='<minimum UMIs per cell>',
                        help='minimum total UMIs required to keep a cell')

    args = parser.parse_args()
    return args

def hammingDist(x, y):
    # find the Hamming distance between two input strings:
    if len(x)!=len(y):
        hd = len(x)
    else:
        hd = 0
        for i in range(len(x)):
            if x[i]!=y[i]:
                hd+=1    # count mismatches
    return hd

def testUmiDist(singList, multList):
    ## test each element of singList against all elements of multList to see if any
    ## singletons are one base off from any of the multis.

    badUmis = []   # output list of UMIs to be removed
    repUmis = {}   # UMIs to be incremented, with the increment count

    for bc0 in singList:
        bc1off = []   # barcodes with a Hamming distance of 1 from the singelton
        for bc1 in multList:
            if hammingDist(bc0, bc1)==1:
                bc1off.append(bc1)  # add to the list of one-offs
        # if a barcode is one off from any multi, remove it:
        if len(bc1off)>0:
            badUmis.append(bc0)
        # if bc0 is one off from only ONE multi, increment the count for that multi:
        if len(bc1off)==1:
            repUmis.setdefault(bc1off[0],0)
            repUmis[bc1off[0]]+=1

    return [badUmis, repUmis]
    
def sumValues(d):
    ## sums the values of the values in a key:value list:
    x=0
    for key,value in d:
        x+=value
    return x

def run(args):
    inFile = args.i      # input file name
    outFile = args.o     # output file name
    readsFile = args.r   # optional output reads file
    nMin = args.n        # minimum number of UMI counts for merging with singletons
    uMin = args.u        # minimum number of UMIs per cell

    print('inFile: %s\n' % inFile)
    print('outFile: %s\n' % outFile)
    print('readsFile: %s\n' % readsFile)

    fIn = open(inFile, 'r')
    umiDict = {}
    umiHist = {}

    gDict = {}   # dictionary whose keys are all observed genes
    nLines = 0

    while 1:
        line = fIn.readline()
        if not line:
            break

        if line.endswith('\n'):
            line=line[:-1]

        nLines+=1
        # split into barcode and count:
        fields = line.split('\t')

        # signal bad input line:
        if len(fields)<4:
            print 'Too few items in line: %s' % line

        # parse the line:
        bc = fields[0]
        g = fields[1]
        umi = fields[2]
        n = int(fields[3])

        # update the gene dictionary:
        gDict.setdefault(g,0)

        # UMI counts histogram:
        umiHist.setdefault(bc,{})
        umiHist[bc].setdefault(g,{})
        umiHist[bc][g].setdefault(umi,0)
        umiHist[bc][g][umi]+=n
            
        # keep stats on all UMIs:
        umiDict.setdefault(umi,0)
        umiDict[umi]+=n

    fIn.close()

    ## for bulk data, we don't do the singlet UMI removal since it's more likely to see a 1-hamming distance UMI that is a true singlets
    ## get lists of UMIs with one count and UMIs with >nMin counts:
    #meanMin = 2 # minimum number of UMI counts for use in computing UMI mean
    #for bc in umiHist.keys():
    #    for g in umiHist[bc].keys():
    #        # for average UMI count calculations:
    #        nzUmis = 0   # number of UMIs with non-zero counts
    #        ntUmis = 0   # number of UMIs with >1 count
    #        umiSum = 0   # sum for UMI mean calculation
    #        # make lists of singlet and multis:
    #        singlets = []
    #        multis = []
    #        for umi in umiHist[bc][g].keys():
    #            if umiHist[bc][g][umi]==1:
    #                singlets.append(umi)
    #            elif umiHist[bc][g][umi]>=nMin:
    #                multis.append(umi)
    #                nzUmis+=1
    #            else:
    #                nzUmis+=1
    #
    #        # separate true singlets from ones that have a Hamming distance of 1 from one or more of the multis:
    #        [rmUmis, incUmis] = testUmiDist(singlets, multis)
    #
    #        # remove bad UMIs from the dictionary
    #        for i in range(len(rmUmis)):
    #            del umiHist[bc][g][rmUmis[i]]     # delete error UMI
    #        # update the counts for UMIs that were uniquely one-off from one of the singletons:
    #        for u in incUmis.keys():
    #            umiHist[bc][g][u]+=incUmis[u]
    #
    ## For each gene, the exprssion value is len(umiHist[g].keys()) (i.e., the number of UMIs 
    ## remaining with at least one copy):
    fOut = open(outFile, 'w')

    ## open the optional read counts file:
    if readsFile==None:
        writeReads=False
    else:
        writeReads=True
        fReads = open(readsFile, 'w')

    ## remove any barcodes/cells with fewer than uMin total UMIs (if uMin>0):
    if uMin>0:
        bcRm = []
        for bc in umiHist.keys():   # loop over barcodes
            bcSum = 0
            for g in umiHist[bc].keys():   # loop over genes
                bcSum += len(umiHist[bc][g].keys())  # sum total UMIs for this cell
            ## remove this barcode if the total count is less than uMin:
            if bcSum<uMin:
                bcRm.append(bc)
    ## remove low-count cells:
    for bc in bcRm:
        del umiHist[bc]
            
    ## file header:
    hStr = 'gene'      # header
    bcList = umiHist.keys()     # just to make absolutely sure they stay in the correct order
    for bc in bcList:
        hStr = '%s\t%s' % (hStr,bc)
    fOut.write('%s\n' % hStr)

    ## option reads count file header
    if writeReads:
        fReads.write('%s\n' % hStr)

    # per-gene count:
    for g in gDict.keys():  # loop over all observed genes
        oStr = '%s' % g
        rStr = '%s' % g
        for bc in bcList:   # UMI counts for each cell
            if umiHist[bc].has_key(g):
                oStr = '%s\t%d' % (oStr,len(umiHist[bc][g].keys()))
                rStr = '%s\t%d' % (rStr,sumValues(umiHist[bc][g].items()))
            else:            
                oStr = '%s\t0' % oStr
                rStr = '%s\t0' % rStr

        fOut.write('%s\n' % oStr)        
        if writeReads:
            fReads.write('%s\n' % rStr)        

    fOut.close()
    if writeReads:
        fReads.close()

    return                      

if __name__ == "__main__":
    args = get_args()
    run(args)


