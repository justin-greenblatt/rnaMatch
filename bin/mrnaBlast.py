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



class Parser:
    def __init__(self, mrnaDirectorie):

        self.mrnaDirectorie = mrnaDirectorie
        self.mrnaParser = SeqIO.parse(mrnaDirectorie, "fasta")
        self.genomeId = mrnaDirectorie.split('/')[-1].rstrip(".gff")
        self.dirDict = {
            "mrnaFolder" : join(dConfig["resources"]["mrna_folder"], self.genomeId),
            "blastDBFolder" : join(dConfig["resources"]["mrna_blast_db_folder"], self.genomeId),
            "blastTestOutDir" : join(dConfig["resources"]["mrna_blast_test_out_folder"], self.genomeId),
            "blastControlOutDir" : join(dConfig["resources"]["mrna_blast_control_out_folder"], self.genomeId)
              }


    def getMrna(self):
        return next(self.mrnaParser, None)

class BlastMrna:

    """
    A class to initiate a revBlast process and hold its directories. Contructer uses a gtf gene entrie.
    from that it opens handlers to output a specific genomewalk and genomeWalkControl csv file. It also
    extracts the seqence of the file form a GenomeDict instance and saves is in a new fasta file. The 
    last function of the class in getRevBlastCommand() and getRevBlastControlCommand() builds a command
    for Popen to open a revBlast.py process.
    """

    #Defining variables shared by instances. Folders and GenomeDict instance.

    def __init__(self, mrna, dirDict):
        self.dirDict = dirDict
        self.mrnaId = mrna.id
        self.fastaFileDir = join(self.dirDict["mrnaFolder"], f"{self.mrnaId}.fa")

        if os.path.isfile(self.fastaFileDir):
            print(f"{self.mrnaId} already exists")
            copies = len(list([a for a in os.listdir(dConfig["resources"]["mrna_folder"]) 
                              if a.split('.')[0].startswith(self.mrnaId)]))
            self.mrnaId = "{}({})".format(self.mrnaId, copies)
            self.fastaFileDir = join(self.dirDict["mrnaFolder"], f"{self.mrnaId}.fa")

        #Defining filenames and directories for species resources.
        self.blastDBDir = join(self.dirDict["blastDBFolder"], f"{self.mrnaId}")
        self.blastTestOutDir = join(self.dirDict["blastTestOutDir"], f"{self.mrnaId}.xml")
        self.blastControlOutDir =  join(self.dirDict["blastControlOutDir"], f"{self.mrnaId}.xml")
        fastaHandler = open(self.fastaFileDir, 'w')
        fastaHeader = ">" + self.mrnaId + "\n"
        fastaHandler.write(fastaHeader)
        fastaHandler.write(str(mrna.seq) + '\n')
        fastaHandler.close()

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

    def run(self):

        self.createBlastDBProcess = Popen(self.getBlastDBCommand(), stdin = PIPE, stdout = PIPE)
        self.blastnStarted = False

    def poll(self) -> int:
        if self.createBlastDBProcess.poll() == 0:
            if self.blastnStarted == False:
                self.revBlastProcess = Popen(self.getTestCommand(), stdin = PIPE, stdout = PIPE)
                self.revBlastControlProcess = Popen(self.getControlCommand(), stdin = PIPE, stdout = PIPE)
                self.blastnStarted = True
                return 1
            else: 
                a = self.revBlastProcess.poll()
                b = self.revBlastControlProcess.poll()
                if a == 0 and b == 0:
                    return 0
                else:
                    return 1
        else:
            return 1

    def close(self):
        #needs to be done to close stdin or stdout files so overall number of open files does not overflow
        do,de = self.createBlastDBProcess.communicate()
        to,te = self.revBlastProcess.communicate()
        co,ce = self.revBlastControlProcess.communicate()
 
class BlastProcessStack:

    """
    A class for managing paralel processes of revBlast.
    """
    def __init__(self):
        self.stackSize = int(pConfig["genomeWalk"]["MAX_PROCESSES"])
        self.refreshTime = float(pConfig["genomeWalk"]["REFRESH_STACK_SEC"])
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
                    print(f"{time.time()} removing {p.mrnaId} , adding {newP.mrnaId}")
                    self.processes.append(newP)
            time.sleep(self.refreshTime)
        print("stack empty. finishing last processes")
        while any(list([a.poll() != 0 for a in self.processes])):
            time.sleep(self.refreshTime)
        self.finished.extend(self.processes)
        self.processes = []

class MrnaBlastExperiment:
    def __init__(self, mrnaDir : str):
        
        self.stack = BlastProcessStack()
        self.mrnaParser = SeqIO.parse(mrnaDir, "fasta")
        self.genomeId = mrnaDir.split('/')[-1].rstrip(".fna")
        self.dirDict = {
            "mrnaFolder" : join(dConfig["resources"]["mrna_genes_folder"], self.genomeId),
            "blastDBFolder" : join(dConfig["resources"]["mrna_blast_db_folder"], self.genomeId),
            "blastTestOutDir" : join(dConfig["resources"]["mrna_blast_test_out_folder"], self.genomeId),
            "blastControlOutDir" : join(dConfig["resources"]["mrna_blast_control_out_folder"], self.genomeId)
              }
        for d in self.dirDict.values():
            if not os.path.isdir(d):
                os.mkdir(d)

    def runExperiment(self):

        print("initializing stack")
        
        geneCount = 0
        for s in self.mrnaParser:
            
            self.stack.add(BlastMrna(s, self.dirDict))
            geneCount += 1
            
        print(f"starting paralel processing of {geneCount} genes")
        self.stack.execute()
