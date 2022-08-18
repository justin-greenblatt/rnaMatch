
"""
Developed for python3.9
justingreeblatt@github | last updated 3/08/2022
dependacies are :
python libraries:
biopython==1.78

CLI:
ncbiblast+
"""
import time
from subprocess import Popen, PIPE
from settings import pConfig, dConfig, loggingConfigPath
import os
import sys
from Bio import SeqIO
import re
import gzip
import logging
import logging.config
from os.path import join
#Setting up Logging
logging.config.fileConfig(loggingConfigPath)
logger = logging.getLogger(__name__)

class GenomeDict:
    """
    A class for accessing a fasta genome. It envelops the biopython SeqIO and allows the program
    to access a specific sequence of a chromossome through its getSeq() method. Its constructer
    only needs a valid fasta genome directorie. At any time only one chromossome sequence is 
    hold in memory.
    """
    def __init__(self, genomeDirectorie):

        self.genomeDirectorie = genomeDirectorie
        self.genomeDict = SeqIO.index(genomeDirectorie, "fasta")
        self.keys = list(self.genomeDict.keys())
        self.chrom = self.genomeDict[self.keys[0]]
     

    def extractSeq(self, geneStart : int, geneEnd : int) -> str:
        """
        A helper function for extracting sequence out of the current chromossom hold in memory
        """
        return str(self.chrom.seq[max(0, geneStart - int(pConfig["genomeWalk"]["DOWNSTREAM_GENE_FLANK"])):min(geneEnd + int(pConfig["genomeWalk"]["UPSTREAM_GENE_FLANK"]), len(self.chrom.seq))])



    def getSeq(self, chromId : str, geneStart : int, geneEnd: int) -> str:
        """
        Get a sequence for a specified chromossome while minimizing loading new chromossomes to memory
        """
        seqStart = min((geneStart, geneEnd))
        seqEnd = max((geneStart, geneEnd))
        if chromId == self.chrom.id:
            return self.extractSeq(seqStart, seqEnd)
        else:
            if chromId in self.keys:
                self.chrom = self.genomeDict[chromId]
                print(f"switching to chromossome {chromId}")
                return self.extractSeq(seqStart, seqEnd)
            else:
                print(f"chromossome {chromId} not in {self.genomeDirectorie}. skipping entrie") 
                return None


class revBlastGene:

    """
    A class to initiate a revBlast process and hold its directories. Contructer uses a gtf gene entrie.
    from that it opens handlers to output a specific genomewalk and genomeWalkControl csv file. It also
    extracts the seqence of the file form a GenomeDict instance and saves is in a new fasta file. The 
    last function of the class in getRevBlastCommand() and getRevBlastControlCommand() builds a command
    for Popen to open a revBlast.py process.
    """

    #Defining variables shared by instances. Folders and GenomeDict instance.
    genomeId = IN_FILE_GENOME.split('/')[-1].rstrip(".fna")
    genomeDict = GenomeDict(IN_FILE_GENOME)
    dirDict = {
            "preMrnaFolder" : join(dConfig["resources"]["pre_mrna_folder"], genomeId),
            "genomeWalkTestFolder" : join(dConfig["resources"]["genome_walk_folder"], genomeId),
            "genomeWalkControlFolder" : join(dConfig["resources"]["genome_walk_control_folder"], genomeId),
            "blastDBFolder" : join(dConfig["resources"]["blast_db_folder"], genomeId),
            "blastTestOutDir" : join(dConfig["resources"]["blast_test_out_folder"], genomeId),
            "blastControlOutDir" : join(dConfig["resources"]["blast_control_out_folder"], genomeId)
              }

    def __init__(self, gtfRecord, genome):
        self.genome = genome
        self.gtfRecord = gtfRecord
        # checking if genome/species specific folders exist.
        for d in self.__class__.dirDict.values():
            if not os.path.isdir(d):
                try:
                    os.mkdir(d)
                except:
                    logging.error("File system corrupted or no permissions to create dir")
                    sys.exit(1)
        #extracting data from gtf entrie
        self.geneId = re.search(r'gene_id \"(.+?)\"', gtfRecord).group(1)
        self.chromId = gtfRecord.split('\t')[0]
        self.geneStart = int(gtfRecord.split('\t')[3])
        self.geneEnd = int(gtfRecord.split('\t')[4])
        self.geneStrand = gtfRecord.split('\t')[6]
        #print(f"initiating {self.geneId}") 
        #Define fasta file for extracted gene. If gene already exists. resolve name clash.
        self.fastaFileDir = join(self.__class__.dirDict["preMrnaFolder"], f"{self.geneId}.fa")
        if os.path.isfile(self.fastaFileDir):
            print(f"{self.geneId} already exists")
            copies = len(list([a for a in os.listdir(dConfig["resources"]["pre_mrna_folder"]) 
                              if a.split('.')[0].startswith(geneId)]))
            self.geneId = "{}({})".format(self.geneId, copies)
            self.fastaFileDir = join(self.__class__.dirDict["preMrnaFolder"], f"{self.geneId}.fa")
        #Defining filenames and directories for species resources.
        self.genomeWalkFileDir = join(self.__class__.dirDict["genomeWalkTestFolder"], f"{self.geneId}.csv")
        self.genomeWalkControlFileDir = join(self.__class__.dirDict["genomeWalkControlFolder"], f"{self.geneId}_control.csv")
        self.blastDBDir = join(self.__class__.dirDict["blastDBFolder"], f"{self.geneId}")
        self.blastTestOutDir = join(self.__class__.dirDict["blastTestOutDir"], f"{self.geneId}")
        self.blastControlOutDir =  join(self.__class__.dirDict["blastControlOutDir"], f"{self.geneId}")

    def getBlastDBCommand(self) -> list[str]:

        command = list([a for a in 
                    [pConfig["makeblastdb"]["RUN_PATH"],
                    "-in" , self.fastaFileDir, 
                    "-out", self.blastDBDir,
                    "-dbtype", pConfig["makeblastdb"]["DB_TYPE"],
                    pConfig["makeblastdb"]["EXTRA_FLAGS"],
                    pConfig["makeblastdb"]["EXTRA_FLAG_VALUE_PAIRS"]]
                    if a])

        return command

    def getTestCommand(self):

        command = list([a for a in 
                    [pConfig["blastn"]["RUN_PATH"],
                    "-query", self.fastaFileDir, 
                    "-strand", "minus",
                    "-db", self.blastDBDir,
                    "-out", self.blastTestOutDir,
                    "-outfmt", pConfig["blastn"]["OUT_FORMAT"],
                    pConfig["blastn"]["EXTRA_FLAGS"],
                    pConfig["blastn"]["EXTRA_FLAG_VALUE_PAIRS"]]
                    if a])

        return command

    def getControlCommand(self) -> list[str]:
    
        command = list([a for a in 
                    [pConfig["blastn"]["RUN_PATH"],
                    "-query", self.fastaFileDir, 
                    "-strand", "plus",
                    "-db", self.blastDBDir,
                    "-out", self.blastControlOutDir,
                    "-outfmt", pConfig["blastn"]["OUT_FORMAT"],
                    pConfig["blastn"]["EXTRA_FLAGS"],
                    pConfig["blastn"]["EXTRA_FLAG_VALUE_PAIRS"]]
                    if a])

        
        return command

    def createFiles(self):
 
        #write out fasta file of extracted sequence.
        fastaHandler = open(self.fastaFileDir, 'w')
        fastaHeader = ">" + ",".join([self.geneId, self.chromId, str(self.geneStart), str(self.geneEnd), str(self.geneStrand)]) + "\n"
        fastaHandler.write(fastaHeader)
        seq = self.__class__.genomeDict.getSeq(self.chromId, self.geneStart, self.geneEnd)
        fastaHandler.write(seq + '\n')
        fastaHandler.close()

        #write out headers for genomeWalk csv files.
        gwHandler = open(self.genomeWalkFileDir, 'w')
        gwcHandler = open(self.genomeWalkControlFileDir, 'w')
        gwHeader = "geneId,chrom,geneStart,geneEnd,geneStrand,queryStart,queryEnd,subjectStart,subjectEnd,matchlength,matchPct\n"
        gwHandler.write(gwHeader)
        gwcHandler.write(gwHeader)
        gwcHandler.close()
        gwHandler.close()

    def run(self):
        self.createFiles()
        self.createBlastDBProcess = Popen(self.getBlastDBCommand(), stdin = PIPE, stdout = PIPE)
        self.createBlastDBProcess.wait()
        self.revBlastProcess = Popen(self.getTestCommand(), stdin = PIPE, stdout = PIPE)
        self.revBlastControlProcess = Popen(self.getControlCommand(), stdin = PIPE, stdout = PIPE)

    def poll(self) -> int:
        a = self.revBlastProcess.poll()
        b = self.revBlastControlProcess.poll()
        if a == 0 and b == 0:
            return 0
        else:
            return 1
    def close(self):
        #needs to be done to close stdin or stdout files so overall number of open files does not overflow
        do,de = self.createBlastDBProcess.communicate()
        to,te = self.revBlastProcess.communicate()
        co,ce = self.revBlastControlProcess.communicate()
 
class RevBlastProcessStack:

    """
    A class for managing paralel processes of revBlast.
    """
    def __init__(self):
        self.stackSize = int(pConfig["genomeWalk"]["MAX_PROCESSES"])
        self.refreshTime = int(pConfig["genomeWalk"]["REFRESH_STACK_SEC"])
        self.stack = []
        self.processes = []
        self.finished = []

    def add(self, revBlastObj):
        self.stack.append(revBlastObj)

    def execute(self):
        print("fulling up process stack for first time")
        for l in range(self.stackSize):
            if len(self.stack) != 0:
                newP = self.stack.pop(0)
                newP.run()
                self.processes.append(newP)
        print("cycling fineshed processes. moving finished to self.finished and adding new from self.stack")
        while len(self.stack) >= len(self.processes):
            for i,p in enumerate(self.processes):
                if p.poll() == 0:
                    oldP = self.processes.pop(i)
                    oldP.close()
                    self.finished.append(oldP)
                    newP = self.stack.pop(0)
                    newP.run()
                    print(f"{time.time()} removing {p.geneId} , adding {newP.geneId}")
                    self.processes.append(newP)
            time.sleep(self.refreshTime)
        print("stack empty. finishing last processes")
        while any(list([a.poll() != 0 for a in self.processes])):
            time.sleep(self.refreshTime)
        self.finished.extend(self.processes)
        self.processes = []

def main():
    #initiating stack
    print("initializing stack")
    rbStack = RevBlastProcessStack()
    gtfHandler = open(IN_FILE_GTF,'r')
    geneCount = 0
    for line in gtfHandler:
        if line.startswith('#') or line.split('\t')[2] not in pConfig["genomeWalk"]["GTF_KEYS"].split(','):
            pass
        else:
            rbStack.add(revBlastGene(line))
            geneCount += 1
            #print(geneCount)
    print(f"starting paralel processing of {geneCount} genes")
    rbStack.execute()

if __name__ == "__main__":
    main()
