
"""
Developed for python3.9
justingreeblatt@github | last updated 3/08/2022
dependacies are :
python libraries:
biopython==1.78

CLI:
ncbiblast+
"""

from subprocess import Popen, PIPE
from settings import pConfig, dConfig, loggingConfigPath
import os
import sys
from Bio import SeqIO
import re
import gzip
import logging
import logging.config

#Setting up Logging
logging.config.fileConfig(loggingConfigPath)
logger = logging.getLogger(__name__)

#Genome
IN_FILE_GENOME = sys.argv[1]

#GTF
IN_FILE_GTF = sys.argv[2]

logger.debug(f"STARTING genomeWalk:{sys.argv[0]} {sys.argv[1]} {sys.argv[2]} {sys.argv[3]} {sys.argv[4]}")


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
        self.keys = list(genomeDict.keys())
        self.chrom = None
     

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
        if self.chrom != None:
            if chromId == self.chrom.id:
                return self.extractSeq(seqStart, seqEnd)
            else:
                if chromId in self.keys:
                    self.chrom = self.genomeDict[chromId]
                    logging.debug(f"switching to chromossome {chromId}")
                    return self.extractSeq(seqStart, seqEnd)
                else:
                    logging.warning(f"chromossome {chromId} not in {self.genomeDirectorie}. skipping entrie") 
                    return None
        elif chromId in self.keys:
            self.chrom = self.genomeDict[chromId]
            logging.debug(f"switching to chromossome {chromId}")
            return self.extractSeq(seqStart, seqEnd)
        else:
            logging.warning(f"chromossome {chromId} not in {self.genomeDirectorie}. skipping entrie") 
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
    genomeId = IN_FILE_GENOME.split('_')[0])
    preMrnaFolder = join(dConfig["resources"]["pre_mrna_folder"], genomeId)
    genomeWalkTestFolder = join(dConfig["resources"]["genomewalk_folder"], genomeId)
    genomeWalkcontrolFolder = join(dConfig["resources"]["genomewalk_control_folder"], genomeId)
    genomeDict = GenomeDict(IN_FILE_GENOME)

    def __init__(self, gtfRecord):

        self.gtfRecord = gtfRecord
        # checking if genome/species specific folders exist.
        for d in [self.__class__.preMrnaFolder,
                  self.__class__.genomeWalkTestFolder,
                  self.__class__.genomeWalkControlFolder]:
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
       
        #Define fasta file for extracted gene. If gene already exists. resolve name clash.
        self.fastaFileDir = join(self.__class__.preMrnaFolder, f"{self.geneId}.fa")
        if os.path.isfile(self.fastaFileDir):
            logging.warning(f"{self.geneId} already exists")
            copies = len(list([a for a in os.listdir(dConfig["resources"]["pre_mrna_folder"] 
                              if a.split('.')[0].startswith(geneId)]))
            self.geneId = "{}({})".format(self.geneId, copies)
            self.fastaFileDir = join(self.__class__.preMrnaFolder, f"{self.geneId}.fa")

        self.genomeWalkFileDir = join(self.__class__.genomeWalkTestFolder, f"{self.geneId}.csv")
        self.genomeWalkControlFileDir = join(self.__class__.genomeWalkControlFolder, f"{self.geneId}_control.csv")

        #write out fasta file of extracted sequence.
        fastaHandler = open(self.fastaFileDir, 'w')
        fastaHeader = ">" + ",".join([self.geneID, self.chromId, str(self.geneStart), str(self.geneEnd), str(self.geneStrand)]) + "\n"
        fastaHandler.write(fastaHeader)
        seq = self.__class__.genomeDict.getSeq(self.chromId, self.geneStart, self.geneEnd)
        fastaHandler.write(seq + '\n')
        fastaHandler.close()
        #write out headers for genomeWalk csv files.
        gwHandler = open(self.genomeWalkFileDir, 'w')
        gwcHandler = open(self.genomeWalkFileControl, 'w')
        gwHeader = "geneId,chrom,geneStart,geneEnd,geneStrand,queryStart,queryEnd,subjectStart,subjectEnd,matchlength,matchPct\n"
        gwHandler.write(gwHeader)
        gwcHandler.write(gwHeader)


    def getTestCommand(self) -> List[str]:
    """
    function for generating a revBlast.py test command for Popen to open a process later.
    """
        preCommand  = list([a for a in [
                pConfig["revBlast"]["USER"],
                pConfig["revBlast"]["USER_FLAGS"],
                pConfig["revBlast"]["INTERPRETER"],
                pConfig["revBlast"]["INTERPRETER_FLAGS"]]
                        if a])

        testCommand = preCommand + [dConfig["scripts"]["REV_BLAST_PATH"], self.fastaFileDir, self.genomeWalkFileDir, "minus"]
        return testCommand

    def getControlCommand(self) -> List[str]:
        """
        function for generating a revBlast.py control commanf for Popen to open a process later.
        """
        preCommand  = list([a for a in [
                pConfig["revBlast"]["USER"],
                pConfig["revBlast"]["USER_FLAGS"],
                pConfig["revBlast"]["INTERPRETER"],
                pConfig["revBlast"]["INTERPRETER_FLAGS"]]
                        if a])

        controlCommand = command + [dConfig["scripts"]["REV_BLAST_PATH"], self.fastaFileDir, self.genomeWalkControlFileDir, "plus"]
        return controlCommand

    def run(self):
        self.revBlastProcess = Popen(self.getTestCommand(), stdin = PIPE, stdout = PIPE)
        self.revBlastControlProcess = Popen(self.getControlCommand(), stdin = PIPE, stdout = PIPE)

    def poll(self):
        a = self.revBlastProcess.poll()
        b = self.revBlastProcess.poll()
        if a == 0 and b == 0:
            return 0
        else:
            return 1

class RevBlastProcessStack:
    """
    A class for managing paralel processes of revBlast.
    """
    def __init__(self):
        self.stackSize = pConfig["genomeWalk"]["MAX_PROCESSES"]
        self.refreshTime = int(pConfig["genomeWalk"]["REFRESH_STACK_SEC"])
        self.stack = []
        self.processes = []
        self.finished = []

    def add(self, revBlastObj):
        self.stack.append(revBlastObj)

    def execute(self):
        #fulling up process stack for first time
        for l in range(self.stackSize):
            if len(self.stack) != 0:
                newP = self.stack.pop()
                newP.run()
                self.processes.append(newP)
        #cycling fineshed processes. moving finished to self.finished and adding new from self.stack
        while len(self.stack) != 0:
            for i,p in enumerate(self.processes):
                if p.poll() == 0:
                    self.finished.append(self.processes.pop(i))
                    newP = self.stack.pop()
                    newP.run()
                    self.processes.append(newP)
            time.sleap(self.refreshTime)
        while any(list([a.poll() != 0 for a in self.processes])):
            time.sleap(self.refreshTime)
        self.finished.extend(self.processes)
        self.processes = []

    def main():
        rbStack = RevBlastProcessStack()
        gtfHandler = open(gtfHandler,'r')
        geneCount = 0
        for line in gtfHandler:
            if line.startswith('#') or line.split('\t')[2] not in pConfig["genomeWalk"]["GTF_KEYS"].split(','):
                pass
            else:
                rbStack.add(revBlastGene(line))
                geneCount += 1
        logging(f"starting paralel processing of {geneCount} genes")
        rbStack.execute()
