"""
Developed on python3.8
justingreenblatt@github.com | last updated 09/05/2022
The class has the functionality of retrieving NCBI genomic data, holding paths to that data, generating reports, histograms and summaries of the data.
"""

#Standard Library Class and function imports
import gzip
from subprocess import Popen, PIPE
from collections import Counter
from os import chdir, getcwd, remove, listdir
from os.path import join, isfile
from requests import get
from re import findall
from random import choice
from sys import exit

#3rd party class function imports
from Bio import SeqIO
from numpy import histogram, histogram2d, arange

#Standard whole library imports
import json
import logging

#My code imports
from myUtils import downloadFromURL
from Histogram import Histogram, Histogram2d

#My constants/parameters imports
from settings.directories import MRNA_FOLDER, GENOME_FOLDER, GENOME_WALK_FOLDER, RNA_WALK_FOLDER, RNA_WALK_CONTROL_FOLDER, GENOME_WALK_CONTROL_FOLDER, GENOME_WALK_PATH, RNA_WALK_PATH, GTF_FOLDER, REPEAT_MASK_FOLDER, HISTOGRAMS_FOLDER, FILE_DIRECTORIES, LOG_FOLDER, HISTOGRAMS2D_FOLDER, IMAGE_FOLDER
from settings.logs import NCBI_GENOME_LOG_PATH, NCBI_GENOME_LOG_LEVEL


class NcbiData:

    def __init__(self, assembly, linkDict):

        """
        Look for any resources from this assembly present in storage.
        """
        self.species, self.assembly = assembly.split('/')
        self.links = linkDict
        self.fileDirectories = {}
        self.id = assembly.replace("/","_")
        
        #find All resources. The lambda funcion is a dummy function to trick the wrapper into doing its work
        self.updateResources(lambda : 1 -1)
        logging.debug("Creating NCBI data object {}".format(self.id))

    def log(func, logLevel = NCBI_GENOME_LOG_LEVEL):
        #Function for setting the log file and level when running a method of this class
        def logWrapper(self, *args, **kwargs):
            logging.basicConfig(filename=join(NCBI_GENOME_LOG_PATH,self.id + ".log"), level = logLevel)
            func(*args, **kwargs)
        
        return  logWrapper
   
    def updateResources(func):
        """
        Decorator for updating resources of this ncbiDataset before and after running an important function
        """
        #Function for looking up for a resource
        def getFile(self,key):
            for f in listdir(RESOURCE_FOLDERS[key]):
                if self.assembly in f:
                    self.fileDirectories[key] = f
                    logging.debug("Found resource {} : {}".format(self.fileDirectories[key], key))
        #find all resource described in RESOURCE_FOLDERS
        def findResources(self):
            for k in RESOURCE_FOLDERS:
                getFile(self,k)

        #decorator/wrapper magic
        def wrapper(self, *args, **kwargs):

            findResources(self)
            logging.debug("Refreshing resources before running object function - {}".format(self.id))
            func(self, *args, **kwargs)
            logging.debug("Refreshing resources after running object function - {}".format(self.id))
            findResources(self)

        return wrapper


    @log
    @updateResources
    def getResource(self, resourceName, resourceRegex, resourceFolder):
        """
        This is a general function for downloading a file from the ncbi Ftp service
        """
        logging.debug("Getting resource {} for {} and storing at {}".format(resourceName, self.id, resourceFolder))
        #Check if resource alreadu exists.
        if not resourceName in self.fileDirectories:

            #Get resource link from links
            resourceLink = [b for a,b in self.links if re.match(resourceRegex, a)][0]
            
            #save current location
            old = getcwd()
            
            #go to Resource folder
            chdir(RESOURCE_FOLDERS[resourceName])
            t0 = time()
            #Download file and if compressed gzip format, then decompress. Other compressions not supported
            if resourceLink.endswith(".gz"):
                logging.debug("Downloading and decompressing {}".format(resourceLink,))
                decompressedFile =  join(RESOURCE_FOLDERS[resourceName] ,downloadFromURL(resourceLink), decompress = True)
                self.fileDirectories[resourceName] = decompressedFile
                logging = logging.debug("Downloaded and Decompressed File : {} ; deltaT : {}".format(decompressedFile, time() -t0))
            else: 
                logging.debug("Downloading {}".format(resourceLink))
                self.fileDirectories[resourceName] = join(RESOURCE_FOLDERS[resourceName] ,downloadFromURL(resourceLink))
                logging.debug("Downloaded  {}".format(resourceLink, time()))
            #go back to old folder
            chdir(old)
        else:
            logging.debug("Attempted to download {} but file already exists at {}".format(resourceName, self.fileDirectories[resourceName]))
    
    def getGenome(self):
        #Look for fasta genome associated to this object. If does not exist. Download genome from ensembl FTP

        self.getResource("dna", r".*(?<!from)_genomic.fna.gz$")

    def getGtf(self):
        
        #Look for fasta genome associated to this object. If does not exist. Download genome from ensembl FTP
          
        self.getResource("gtfAnnotation", r".*(?<!from)_genomic.gtf.gz$")

    def getMrna(self):

        #Look for fasta rna associated to this object. If does not exist. Download rna from ncbi FTP
          
        self.getResource("mrna", r".*rna_from_genomic.gtf.gz$")

    def generateHistograms(self):
        """Generate histograms of data associated to this object.
        """

        logging.info("generating histograms for {}|{}".format(self.species, self.assembly))

        self.histograms = {}
        self.histograms2d = {}
        
        if "genomeWalk" in self.fileDirectories:
        
            self.histograms["gwSizeHist"] = Histogram(self, "genomeWalk", "gwSizeHist",
                lambda x: int(x.split(',')[-2]),
                1000, (float(0),float(1000)))
        
            self.histograms["gwCSizeHist"] = Histogram(self,"genomeWalkControl", "gwCSizeHist",
                          lambda x: int(x.split(',')[-2]),
                          1000, (0,1000))
        
        
            self.histograms2d["gwSizeHist2d"] = Histogram2d(self, "genomeWalk",
                          "gwSizeHist2d",
                          lambda x: list([float(a) for a in x.split(',')[-2:]]),
                          (500,35),
                          ((0, 500), (0.65, 1)))

            self.histograms2d["gwCSizeHist2d"] = Histogram2d(self, "genomeWalkControl",
                          "gwCSizeHist2d",
                          lambda x: list([float(a) for a in x.split(',')[-2:]]),
                          (500,35),
                          ((0, 500), (0.65, 1)))
        else:
            logging.warning("No GenomeWalkFile for {}|{}".format(self.species, self.assembly))

        #These next 3 lines are a fix for a bug
        if "mSizeHist" in self.fileDirectories:
            logging.info("fixing the mSizeHist bug for {}|{}".format(self.genome.species, self.genome.assembly))
            self.fileDirectories["rmSizeHist"] = self.fileDirectories["mSizeHist"]
            self.fileDirectories.pop("mSizeHist")

        if "repeatMask" in self.fileDirectories:
            self.histograms["rmSizeHist"] = Histogram(self, "repeatMask", "rmSizeHist",
                                              lambda x: (lambda z: int(z[1]) - int(z[0]))([a for a in x.split() if a][5:7]),
                                                                        1000, (0, 1000), 3)
        else:
            logging.warning("No repeatMask for {}|{}".format(self.species, self.assembly))

        logging.debug("generated histograsm for ncbi object: name species histograms histograms2d [{},{},{},{}]".format(self.species, self.assembly, str(self.histograms), str(self.histograms2d)))

    @log
    @updateResources
    def runGenomeWalk(self):
        """Run genomeWalk.py for this object. Download input gtf and genome files if they are not downloaded yet. Then create process
        """
        logging.info("Running genomeWalk for {}|{}".format(self.assembly, self.species))

        blastOutFile = join(GENOME_WALK_FOLDER, self.id + ".gw")
        blastControlOutFile = join(GENOME_WALK_CONTROL_FOLDER, self.id + "_control.gw")

        if not blastOutFile in self.fileDirectories:
            if not "gtfAnnotation" in self.fileDirectories:
                try:
                    self.getGtf()
                except:
                    logging.error("getGtf failed for {}|{}".format(self.assembly, self.species))

            if not "dna" in self.fileDirectories:
                try:
                    self.getGenome()
                except:
                    logging.error("getGenome failed for {}|{}".format(self.assembly, self.species))

            try:
                 command = ["python3", GENOME_WALK_PATH,
                           self.fileDirectories["dna"],
                           self.fileDirectories["gtfAnnotation"],
                           blastOutFile,
                           blastControlOutFile]

                 p = Popen(command,stdout = PIPE)
                 logging.debug("Starting process genomeWalk: command id time [{},{},{}]".format(command, p.id, time()))
                 p.wait()
                 logging.debug("Finished process genomeWalk: command id time [{},{},{}]".format(command, p.id, time()))
                 remove(self.fileDirectories["dna"])
                 remove(self.fileDirectories["gtfAnnotation"])
                 self.fileDirectories.pop("gtfAnnotation")
                 self.fileDirectories.pop("dna")
                 logging.debug("Removed gtf and genome file: genomeFile gtfFile [{},{}]".format(self.fileDirectories["dna"],self.fileDirectories["gtfAnnotation"]))
            except:
                 logging.error("failed to run genomeWalk for {}|{}".format(self.assembly, self.species))

            genomeWalkFiles = listdir(GENOME_WALK_FOLDER)

            for i in genomeWalkFiles:
                if self.id in i.lower().replace(" ","_") and "control" not in i:
                    self.fileDirectories["genomeWalk"] = join(GENOME_WALK_FOLDER, i)

            for i in genomeWalkFiles:
                if self.id in i.lower().replace(" ","_") and "control" in i:
                    self.fileDirectories["genomeWalkControl"] = join(GENOME_WALK_FOLDER, i)
            logging.debug("Added paths to fileDirectories: genomeWalkPath genomeWalkControlPath [{},{}]".format(self.fileDirectories["genomeWalk"], self.fileDirectories["genomeWalkControl"]))

    @log
    @updateResources
    def runRnaWalk(self):
        """Run rnaWalk.py for this object. Download input gtf and rna files if they are not downloaded yet. Then create process
        """
        logging.info("Running rnaWalk for {}|{}".format(self.assembly, self.species))

        blastOutFile = join(RNA_WALK_FOLDER, self.id + ".gw")
        blastControlOutFile = join(RNA_WALK_CONTROL_FOLDER, self.id + "_control.gw")

        if not blastOutFile in self.fileDirectories:

            if not "mrna" in self.fileDirectories:
                try:
                    self.getMrna()
                except:
                    logging.error("getMrna failed for {}|{}".format(self.assembly, self.species))

            try:
                 command = ["python3", RNA_WALK_PATH,
                           self.fileDirectories["mrna"],
                           blastOutFile,
                           blastControlOutFile]

                 p = Popen(command,stdout = PIPE)
                 logging.debug("Starting process rnaWalk: command id time [{},{},{}]".format(command, p.id, time()))
                 p.wait()
                 logging.debug("Finished process rnaWalk: command id time [{},{},{}]".format(command, p.id, time()))
                 remove(self.fileDirectories["mrna"])
                 remove(self.fileDirectories["gtfAnnotation"])
                 self.fileDirectories.pop("gtfAnnotation")
                 self.fileDirectories.pop("mrna")
                 logging.debug("Removed gtf and rna file: rnaFile gtfFile [{},{}]".format(self.fileDirectories["mrna"],self.fileDirectories["gtfAnnotation"]))
            except:
                 logging.error("failed to run rnaWalk for {}|{}".format(self.assembly, self.species))

            rnaWalkFiles = listdir(RNA_WALK_FOLDER)

            for i in rnaWalkFiles:
                if self.id in i.lower().replace(" ","_") and "control" not in i:
                    self.fileDirectories["rnaWalk"] = join(RNA_WALK_FOLDER, i)

            for i in rnaWalkFiles:
                if self.id in i.lower().replace(" ","_") and "control" in i:
                    self.fileDirectories["rnaWalkControl"] = join(RNA_WALK_FOLDER, i)
            logging.debug("Added paths to fileDirectories: rnaWalkPath rnaWalkControlPath [{},{}]".format(self.fileDirectories["rnaWalk"], self.fileDirectories["rnaWalkControl"]))
