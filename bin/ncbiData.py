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
from os.path import join, isfile, isdir
import xml.etree.ElementTree as xmlParser
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
import premrnaBlast
import mrnaBlast
#My constants/parameters importS
from settings.directories import  RESOURCE_FOLDERS
from settings.resourceLinkRegex import RESOURCE_REGEX
from settings import lConfigPath, dConfig, pConfig
#Setting up Logging
logging.config.fileConfig(lConfigPath)
logger = logging.getLogger("ncbiData")

#Decorators
def updateResources(func : Callable) -> Callable:
    """
    Decorator for updating resources of this ncbiDataset before and after running an important function
    """
    #Function for looking up for a resource
    def getFile(slf, key : str) -> None:
        foundFlag = False
        for f in listdir(dConfig["resources"][key]):
            if slf.assembly in f:
                slf.fileDirectories[key] = join(dConfig["resources"][key],f)
                foundFlag = True
        if not foundFlag:
            slf.fileDirectories.pop(key,"none")

    #find all resource described in dConfig["resources"]
    def findResources(slf):
        for k in dConfig["resources"]:
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

def migrate(func : Callable) -> Callable:
    """
    Decorator for moving new files and folders generated to the shared nfs memmory
    """
    #Function for looking up for a resource
    #decorator/wrapper magic
    def wrapper(slf, *args, **kwargs):

        func(slf, *args, **kwargs)

        for k in slf.fileDirectories:
            originDir = slf.fileDirectories[k]
            destDir = join(dConfig["nfs"][k], slf.fileDirectories[k].split('/')[-1])

            if not (isdir(destDir) or isfile(destDir)):
                p = Popen(["cp", originDir, destDir])
                p.wait()
                logger.debug(f"migrating {originDir} ---> {destDir}\n{p.communicate()}\n")

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

    @migrate
    @updateResources
    def getResource(self, resourceName : str) -> None:
        """
        This is a general function for downloading a file from the ncbi Ftp service
        """
        resourceRegex = RESOURCE_REGEX[resourceName]
        outFolder = dConfig["resources"][resourceName]
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
    @migrate
    @updateResources
    def genBlastReport(self, folderKey, outName):

        hspKeys = ['Hsp_num', 'Hsp_bit-score', 'Hsp_score', 'Hsp_evalue', 'Hsp_query-from',
                   'Hsp_query-to', 'Hsp_hit-from', 'Hsp_hit-to', 'Hsp_identity',
                   'Hsp_positive', 'Hsp_gaps', 'Hsp_align-len']


        summaryOut = open(outName, 'w')
        summaryOut.write("gene_id,chromossome,gene_start,gene_end,gene_strand,")
        summaryOut.write(','.join(list([k.lower() for k in hspKeys])) + '\n')

        bFiles = listdir(self.fileDirectories[folderKey])

        for b in bFiles:
            blastResultsPath = join(self.fileDirectories[folderKey], b)
            xmlHandler = open(blastResultsPath)
            xmlData = xmlParser.parse(xmlHandler)
            r = xmlData.getroot()
            geneData = r.find("BlastOutput_query-def").text + ','
            for hsp in r.iter("Hsp"):
                hspData = ','.join(list([hsp.find(k).text for k in hspKeys])) + '\n'
                summaryOut.write(geneData)
                summaryOut.write(hspData)
            xmlHandler.close()
        summaryOut.close()
    @migrate
    @updateResources
    def runPremrnaBlast(self) -> None:
        if not "gtf" in self.fileDirectories:
            self.getResource("gtf")
        if not "genome" in self.fileDirectories:
            self.getResource("genome")

        
        experiment = premrnaBlast.PremrnaBlastExperiment(self.fileDirectories["genome"], self.fileDirectories["gtf"])
        experiment.runExperiment()
        self.deleteResource("genome")
        self.deleteResource("gtf")
        self.genBlastReport("premrna_blast_test_out", join(dConfig["resources"]["premrna_blast_test_out_summary"],self.id + ".csv"))
        self.genBlastReport("premrna_blast_control_out", join(dConfig["resources"]["premrna_blast_control_out_summary"],self.id + ".csv"))
    @migrate 
    @updateResources
    def runMrnaBlast(self) -> None:

        if not "mrna" in self.fileDirectories:
            self.getResource("mrna")

        
        experiment = mrnaBlast.MrnaBlastExperiment(self.fileDirectories["mrna"])
        experiment.runExperiment()
        self.deleteResource("mrna")
        self.genBlastReport("mrna_blast_test_out", join(dConfig["resources"]["mrna_blast_test_out_summary"],self.id + ".csv"))
        self.genBlastReport("mrna_blast_control_out", join(dConfig["resources"]["mrna_blast_control_out_summary"],self.id + ".csv"))
    @migrate
    @updateResources
    def migrate(self):

        for k in self.fileDirectories:
            originDir = self.fileDirectories[k]
            destDir = os.path.join(dConfig["nfs"][k], self.fileDirectories[k].split('/')[-1])

            if not (os.path.isdir(destDir) or os.path.isfile(destDir)):
                print(f"migrating {originDir} to {destDir}")
                p = Popen(["mv", originDir, destDir])
                p.wait()
"""
    def generateHistograms(self):
        #Generate histograms of data associated to this object.

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
"""
 
