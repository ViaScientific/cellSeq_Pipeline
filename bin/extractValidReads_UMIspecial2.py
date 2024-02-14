import sys
import cPickle
import pysam
import scipy.stats
import random
import getopt
import os.path
from Bio.Seq import Seq
import re
import argparse

## NOTE: this may require python 2.7.5...

BC_LEN = 6    # BC1 length
UMI_LEN = 7    # UMI length
#now modified specially for custom UMI NNVNVNV::includes all UMI
MAX_BCMISMATCH = 1   # number of bases that can be error-corrected

# Argument parsing:
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, metavar='<R1 file>',
                        help='input R1 fastq file')
    parser.add_argument('-o', type=str, required=False,default=None, metavar='<outDir>',
                        help='output file directory')
    parser.add_argument('-b', type=str, required=True, metavar='<validBcFile>',
                        help='valid barcode list file')
    args = parser.parse_args()
    return args

def patMatch(x, y, n):
    # default match location:
    loc = None

    # attempt to find pattern 'x' in query sequence 'y', allowing a maximum of 'n' mismatches:
    lx = len(x)
    ly = len(y)

    for i in range(ly-lx):   # loop over each possible start location
        nm = 0    # number of mismatches found
        for j in range(lx):   # loop over each base of the target pattern
            if x[j]!=y[i+j]:
                nm+=1        # mismatch
                if nm>n:
                    break    # too many mismatches
        if nm<=n:
            loc = i   # match location
            break

    return loc

### Barcode error-correction:
def correctBC(bc, bcDict):
    bcNew = None    # default (uncorrectable)
    nMatch = 0
    for k in bcDict.keys():
        if hammingDist(bc, k)<=MAX_BCMISMATCH:
            bcNew=k   # corrected barcode
            nMatch+=1 # count the number of possible corrections
            ##break     # don't bother looking any further

    if nMatch>1:
        bcNew = None  # barcode can't be unambiguously corrected

    return bcNew   # return corrected barcode

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

def parseBarcodeAndUmi(r1, bcDict):
    ## assumes UMI is the first UMI_LEN bases and the barcode is the next BC_LEN bases
    umiValid = 0        # default=false
    bcValid = 0         # default=false
    bcCorrected = 0     # default=false

    # extract the sequences that should be the UMI and BC:
    umi = r1.sequence[:UMI_LEN]
    bc = r1.sequence[UMI_LEN:(UMI_LEN+BC_LEN)]

    # first, check for valid BC and UMI with no Ns:
    if bcDict.has_key(bc):
        bcValid = 1      # valid BC
    else:
        # check for correctable BC:
        bcCorr = correctBC(bc, bcDict)
        if bcCorr!=None:
            bc = bcCorr       # corrected BC
            bcCorrected = 1   # BC was corrected

    # check for Ns in the UMI:
    if 1:#umi.count('N')==0 and umi[2]!='T' and umi[4]!='T' and umi[6]!='T':
        umiValid = 1   # valid UMI

    return [bc, umi, bcValid, bcCorrected, umiValid]

def fastqWrite(f, r, rName):
    # write each field of the read:
    f.write('@%s %s\n' % (rName, r.comment))     # read name and comment
    f.write('%s\n' % r.sequence)                  # the sequence
    f.write('+%s %s\n' % (rName, r.comment))     # read name and comment (filler?)
    f.write('%s\n' % r.quality)                   # the quality string
    return

##### RUN PROGRAM ######

def run(args):
    basename = args.i       # forward reads file with barcode and UMIs
    outdir = args.o   # output directory
    if outdir==None:
        outdir = os.path.dirname(basename)  # if not provided, put the output file in the same directory as the input
    bcFile = args.b   # valid barcode file

    # load the valid barcode dictionary:
    fIn = open(bcFile,'r')
    bcDict = {}
    while 1:
        line = fIn.readline()
        if not line:
            break
        # skip the header line:
        if line.startswith('well'):
            continue
        if line.endswith('\n'):
            line=line[:-1]

        fields = line.split('\t')

        bc = fields[0]
        bcDict.setdefault(bc,0)

    fIn.close()

    # construct the names of the input files:
    f1 = basename
    f2 = basename.replace('.R1','.R2')   # R1 is the UMI/barcode read
    # and output fastq file:
    fields = f2.split('_')
    for i in range(len(fields)):
        if i==0:
            f2out = fields[i]
        elif i==(len(fields)-1):
            f2out = f2out+'_valid_'+fields[i]
        else:
            f2out = f2out+'_'+fields[i]

    f2out = f2out.replace(".gz","")    # just in case...
    f2out =  os.path.join(outdir, os.path.basename(f2out))

    # open the input file:
    fq1 = pysam.FastqFile(f1)
    fq2 = pysam.FastqFile(f2)
    fq2out = open(f2out, 'w')

    eFlag = False     # error flag
    rCount = 0
    bcValid = 0    # valid barcodes
    bcCorrected = 0    # corrected barcodes
    umiValid = 0   # valid UMI

    countMod = 100000

    # loop over all reads:
    while 1:
        try:
            r1 = fq1.next()     # R1 read
            r2 = fq2.next()     # R2 read
            rCount+=1           # read counter
            #if not rCount%countMod:
            #    print 'read %d' % rCount
        except StopIteration:

            break      # last item
        except:
            print 'pysam.Fastqfile iterator error.'
            eFlag = True
            break

        # get the barcode from the R1 read
        [bc, umi, bcV, bcC, umiV] = parseBarcodeAndUmi(r1, bcDict)

        ## update counts:
        bcValid += bcV
        bcCorrected += bcC
        umiValid += umiV

        ## update counts for this barcode:
        if (bcV or bcC) and umiV:
            cleanName = r2.name.replace('/',':')
            rName = '%s:%s:%s' % (cleanName, bc,umi)   # append <barcode>:<UMI> to R2 read name
            # write the read out to the new R2 fastq file:
            fastqWrite(fq2out, r2, rName)

            # update barcode counts:
            bcDict[bc]+=1

    # close the input files:
    fq1.close()
    fq2.close()
    fq2out.close()

    # print counts:
    print 'Total reads: %d' % rCount
    for bc in bcDict.keys():
        count = bcDict[bc]
        print('\t%s\t%d' % (bc,count))
    print('Valid barcodes:\t%d' % bcValid)
    print('Corrected barcodes:\t%d' % bcCorrected)
    print('Valid UMIs:\t%d' % umiValid)

    return


if __name__ == "__main__":

    args = get_args()
    run(args)
