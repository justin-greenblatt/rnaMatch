"""
Developed on python3.9
justingreenblatt@github.com | last updated 06/06/2022
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
from typing import List, Dict, Callable

#3rd party class function imports
from Bio import SeqIO
from numpy import histogram, histogram2d, arange
from copy import deepcopy
from time import time

#Standard whole library imports
import re
import json
import logging
import logging.config

#My code imports
from myUtils import downloadFromURL, dictComparison
from Histogram import Histogram, Histogram2d

#My constants/parameters importS
from settings.directories import  RESOURCE_FOLDERS, LOGGING_CONF, GENOME_WALK_PATH, RNA_WALK_PATH
from settings.resourceLinkRegex import RESOURCE_REGEX

#Setting up Logging
logging.config.fileConfig(LOGGING_CONF)
logger = logging.getLogger("ncbiData")

#Decorators
def updateResources(func : Callable) -> Callable:
    """
    Decorator for updating resources of this ncbiDataset before and after running an important function
    """
    #Function for looking up for a resource
    def getFile(slf, key : str) -> None:
        foundFlag = False
        for f in listdir(RESOURCE_FOLDERS[key]):
            if slf.assembly in f:
                slf.fileDirectories[key] = join(RESOURCE_FOLDERS[key],f)
                foundFlag = True
        if not foundFlag:
            slf.fileDirectories.pop(key,"none")

    #find all resource described in RESOURCE_FOLDERS
    def findResources(slf):
        for k in RESOURCE_FOLDERS:
            getFile(slf,k)

    #decorator/wrapper magic
    def wrapper(slf, *args, **kwargs):

        findResources(slf)
        oldFiles = deepcopy(slf.fileDirectories)
        func(slf, *args, **kwargs)
        findResources(slf)
        newFiles = deepcopy(slf.fileDirectories)
        logger.debug("Species files comparison after calling {} on{}{}\n {}".format(func.__name__, args, kwargs, dictComparison(oldFiles,newFiles)))

    return wrapper


class ncbiData:

    def __init__(self, assembly : str, linkDict : Dict[str,str]):

        """
        Look for any resources from this assembly present in storage.
        """

        self.species, self.assembly = assembly.split('/')[:2]
        self.links = linkDict
        self.fileDirectories = {}
        self.id = assembly.replace("/","_")
        logger.info("Creating ncbiData object for {}".format(assembly))


        #find All resources. The lambda funcion is a dummy function to trick the wrapper into doing its work
        updateResources(lambda : 1 -1)


    @updateResources
    def getResource(self, resourceName : str) -> None:
        """
        This is a general function for downloading a file from the ncbi Ftp service
        """
        resourceRegex = RESOURCE_REGEX[resourceName]
        outFolder = RESOURCE_FOLDERS[resourceName]
        logger.debug("Getting resource {} for {} and storing at {}".format(resourceName, self.id, outFolder))
        #Check if resource alreadu exists.
        if not resourceName in self.fileDirectories:

            #Get resource link from links
            resourceLinks = [b for a,b in self.links.items() if re.match(resourceRegex, a)]
            resourceLink = resourceLinks[0]
            #save current location
            old = getcwd()
 
            #go to Resource folder
            chdir(outFolder)
            t0 = time()
            #Download file and if compressed gzip format, then decompress. Other compressions not supported
            if resourceLink.endswith(".gz"):
                logger.debug("Downloading and decompressing {}".format(resourceLink,))
                decompressedFile =  join(outFolder ,downloadFromURL(resourceLink, decompress = True))
                self.fileDirectories[resourceName] = decompressedFile
                logger.debug("Downloaded and Decompressed File : {} ; deltaT : {}".format(decompressedFile, time() -t0))
            else:
                logging.debug("Downloading {}".format(resourceLink))
                self.fileDirectories[resourceName] = join(outFolder ,downloadFromURL(resourceLink))
                logging.debug("Downloaded  {}".format(resourceLink, time()))
            #go back to old folder
            chdir(old)
        else:
            logger.debug("Attempted to download {} but file already exists at {}".format(resourceName, self.fileDirectories[resourceName]))

    @updateResources
    def deleteResource(self, resourceName):
        remove(self.fileDirectories[resourceName])


    def generateHistograms(self):
        """Generate histograms of data associated to this object.
        """

        logger.info("generating histograms for {}|{}".format(self.species, self.assembly))

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
            logger.info("fixing the mSizeHist bug for {}|{}".format(self.genome.species, self.genome.assembly))
            self.fileDirectories["rmSizeHist"] = self.fileDirectories["mSizeHist"]
            self.fileDirectories.pop("mSizeHist")

        if "repeatMask" in self.fileDirectories:
            self.histograms["rmSizeHist"] = Histogram(self, "repeatMask", "rmSizeHist",
                                              lambda x: (lambda z: int(z[1]) - int(z[0]))([a for a in x.split() if a][5:7]),
                                                                        1000, (0, 1000), 3)
        else:
            logger.warning("No repeatMask for {}|{}".format(self.species, self.assembly))

        logger.debug("generated histograsm for ncbi object: name species histograms histograms2d [{},{},{},{}]".format(self.species, self.assembly, str(self.histograms), str(self.histograms2d)))


    @updateResources
    def runGenomeWalk(self) -> None:
        """Run genomeWalk.py for this object. Download input gtf and genome files if they are not downloaded yet. Then create process
        """
        logger.info("Running genomeWalk for {}|{}".format(self.assembly, self.species))

        blastOutFile = join(RESOURCE_FOLDERS["genomeWalk"], self.id + ".csv")
        blastControlOutFile = join(RESOURCE_FOLDERS["genomeWalkControl"], self.id + "_control.csv")

        if not blastOutFile in self.fileDirectories:
            if not "gtfAnnotation" in self.fileDirectories:
                try:
                    self.getResource("gtf")
                except:
                    logger.error("getGtf failed for {}|{}".format(self.assembly, self.species))

            if not "dna" in self.fileDirectories:
                try:
                    self.getResource("dna")
                except:
                    logger.error("getGenome failed for {}|{}".format(self.assembly, self.species))

            try:
                 command = ["python3", GENOME_WALK_PATH,
                           self.fileDirectories["dna"],
                           self.fileDirectories["gtfAnnotation"],
                           blastOutFile,
                           blastControlOutFile]

                 p = Popen(command,stdout = PIPE)
                 logger.debug("Starting process genomeWalk: command id [{},{}]".format(command, p.id))
                 p.wait()
                 logger.debug("Finished process genomeWalk: command id [{},{}]".format(command, p.id))
                 self.deleteResource("dna")
                 self.deleteResource("gtf")
            except:
                 logger.error("failed to run genomeWalk for {}|{}".format(self.assembly, self.species))


    @updateResources
    def runRnaWalk(self) -> None:
        """
           Run rnaWalk.py for this object. 
           Download input gtf and rna files if they are not downloaded yet.
           Then create process
        """
        logger.info("Running rnaWalk for {}|{}".format(self.assembly, self.species))

        blastOutFile = join(RESOURCE_FOLDERS["rnaWalk"], self.id + "_rnaWalk.csv")
        blastControlOutFile = join(RESOURCE_FOLDERS["rnaWalkControl"], self.id + "_rnaWalkControl.csv")

        if not blastOutFile in self.fileDirectories:

            if not "mrna" in self.fileDirectories:
                self.getResource("mrna")
            try:
                 command = ["python3", RNA_WALK_PATH,
                           self.fileDirectories["mrna"],
                           blastOutFile,
                           blastControlOutFile]

                 p = Popen(command,stdout = PIPE)
                 logger.debug("Starting process rnaWalk: command id [{},{}]".format(command, p.id))
                 p.wait()
                 logger.debug("Finished process rnaWalk: command id [{},{}]".format(command, p.id))
                 self.deleteResource("mrna")
                 self.deleteResource("gtfAnnotation")
                
            except:
                 logger.error("failed to run rnaWalk for {}|{}".format(self.assembly, self.species))
