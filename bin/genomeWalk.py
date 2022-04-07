#This iterates over the genome and annotation and applies blast on the gene region

from subprocess import Popen, PIPE
import os
import sys
from Bio import SeqIO
import re
from settings.runParameters import MINUS_GENE_REGION, PLUS_GENE_REGION
from settings.directories import REV_BLAST_PATH
import gzip

#Genome
IN_FILE_GENOME = sys.argv[1]

#GTF
IN_FILE_GTF = sys.argv[2]

#outFile
OUT_FILE = sys.argv[3]
CONTROL_OUT_FILE = sys.argv[4]
print("START")

new = open(OUT_FILE,'x')
newControl = open(CONTROL_OUT_FILE, 'x')
new.write("geneId,geneStart,geneEnd,geneStrand,queryStart,queryEnd,subjectStart,subjectEnd,matchlength,matchPct\n")
newControl.write("geneId,geneStart,geneEnd,geneStrand,queryStart,queryEnd,subjectStart,subjectEnd,matchlength,matchPct\n")
new.close()
newControl.close()
genomeHandler = gzip.open(IN_FILE_GENOME,'rt')
genomeIterator = SeqIO.parse(genomeHandler, "fasta")
gtfIterator = gzip.open(IN_FILE_GTF, 'rt')

flag = True

chromosome = next(genomeIterator, None)


def getNextGene(iterator):
    line = next(iterator)
    while line.startswith('#') or line.split('\t')[2] != "gene":
        line = next(gtfIterator, "")
        if line == "":
            sys.exit(0)
            break
    return line

gene = getNextGene(gtfIterator)
#Start Looping over genome and annotations

while flag:
    
    if isinstance(chromosome, type(None)):
        print("END")
        sys.exit(0)
        break
    
    #skipComments
    while gene.split('\t')[0] == chromosome.id:
        
        #Getting data of the gene
        geneID = re.search(r'gene_id \"(.+?)\"', gene).group(1)
        geneStart = int(gene.split('\t')[3])
        geneEnd = int(gene.split('\t')[4])
        geneStrand = gene.split('\t')[6]
        geneSeq = str(chromosome.seq[max(0, geneStart - MINUS_GENE_REGION):min(geneEnd + 2000, len(chromosome.seq))])

        #Writing gene data to temporary file 
        tempFilename = "temp_fasta_" + geneID + '.fa'
        tempGeneFasta = open(tempFilename, 'w')
        fastaHeader = ">" + ",".join([geneID, str(geneStart), str(geneEnd), str(geneStrand)]) + "\n"
        tempGeneFasta.write(fastaHeader)
        tempGeneFasta.write(geneSeq)
        tempGeneFasta.close()
    
        #Running revBlast on gene
        
        command = ["sudo", "python3", REV_BLAST_PATH, tempFilename, OUT_FILE, "minus"]
        p = Popen(command, stdout = PIPE)
        p.wait()
        
        #getControl

        controlCommand = ["sudo", "python3", REV_BLAST_PATH, tempFilename, CONTROL_OUT_FILE, "plus"]
        pc = Popen(controlCommand, stdout = PIPE)
        pc.wait()

        #Removing temporary file
        os.remove(tempFilename)
    
        gene = getNextGene(gtfIterator)
    chromosome = next(genomeIterator, None)
gtfHandler.close()
genomeHandler.close()

