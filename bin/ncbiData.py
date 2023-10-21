"""
Developed on python3.9
justingreenblatt@github.com | last updated 06/06/2022
The class has the functionality of retrieving NCBI genomic data, holding paths to that data, generating reports, histograms and summaries of the data.
"""

#Standard Library Class and function imports
import gzip
import time
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
from shutil import rmtree
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

def migrate(func : Callable, Compressed = True) -> Callable:

    #Decorator for moving new files and folders generated to the shared nfs memmory
    #Function for looking up for a resource
    #decorator/wrapper magic
    def wrapper(slf, *args, **kwargs):

        func(slf, *args, **kwargs)
        print("--------------MIGRATING--------------")
        for k in slf.fileDirectories:
            originDir = slf.fileDirectories[k]
            if originDir.endswith(".gz"):
                destDir = join(dConfig["cloud"][k], originDir.split('/')[-1])

                if (not (isdir(destDir) or isfile(destDir))):
                    p = Popen(["gsutil", "cp", "-r", originDir, destDir])
                    p.wait()
                    logger.debug(f"migrating {originDir} ---> {destDir}\n{p.communicate()}\n")

    return wrapper


class ncbiData:

    def __init__(self, linksDir):

        """
        Look for any resources from this assembly present in storage.
        """
        linksFile = linksDir.split('/')[-1]
        self.species = linksFile.split('-')[0]
        self.id = linksFile.rstrip(".json")
        self.assembly = self.id.lstrip(self.species + '-')
        fileHandler = open(linksDir)
        self.links = json.load(fileHandler)
        fileHandler.close()
        self.fileDirectories = {}
        logger.info("Creating ncbiData object for {}".format(self.assembly))


        #find All resources. The lambda funcion is a dummy function to trick the wrapper into doing its work
        updateResources(lambda : 1 -1)


    @updateResources
    #@migrate
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
                if not isfile(decompressedFile):
                    if isdir(decompressedFile):
                        raise Exception(f"{decompressedFile} is a folder not a file")
                    raise Exception(f"{decompressedFile} does not exist")
            else:
                logging.debug("Downloading {}".format(resourceLink))
                newFile = join(outFolder ,downloadFromURL(resourceLink))
                self.fileDirectories[newFile]
                logging.debug("Downloaded  {}".format(resourceLink, time()))
                if not isfile(newFile):
                    if isdir(newFile):
                        raise Exception(f"{newFile} is a folder not a file")
                    raise Exception(f"{newFile} does not exist")
            #go back to old folder
            chdir(old)
        else:
            logger.debug("Attempted to download {} but file already exists at {}".format(resourceName, self.fileDirectories[resourceName]))

    @updateResources
    def deleteResource(self, resourceName):
        remove(self.fileDirectories[resourceName])

    @updateResources
    @migrate
    def upload(self):
        return

    @updateResources
    def compressAndUpload(self,compressList):

         pool = []
         #uploadList = []
         for name in compressList:
             indir = self.fileDirectories[name]

             if isdir(indir):
                 outdir = indir + ".tar.gz"
                 cprocess = Popen(["tar","-czvf", outdir, indir])
                 #uploadList.append((name,outdir))
                 #pool.append(cprocess)
                 cprocess.wait()
                 rmtree(indir)

             elif isfile(indir):
                 outdir = indir + ".gz"
                 cprocess = Popen([f"bgzip -c {indir} > {outdir} "], shell = True)
                 #uploadList.append((name,outdir))
                 #pool.append(cprocess)
                 cprocess.wait()
                 remove(indir)
             else:
                 raise Exception(f" {name} type file or folder not found, dont forget to implement updateResources to map output to application")
         #while any(list([a.poll() != 0 for a in pool])):
             #time.sleep(2)

             destDir = join(dConfig["cloud"][name], outdir.split('/')[-1])
             p = Popen(["gsutil", "cp", "-r", outdir, destDir])
             p.wait()
    #@migrate
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

    #@migrate
    @updateResources
    def runPremrnaBlast(self) -> None:
        if not "gtf" in self.fileDirectories:
            self.getResource("gtf")
        if not "genome" in self.fileDirectories:
            self.getResource("genome")

        
        experiment = premrnaBlast.PremrnaBlastExperiment(self.fileDirectories["genome"], self.fileDirectories["gtf"])
        experiment.runExperiment()
        self.genBlastReport("premrna_blast_test_out", join(dConfig["resources"]["premrna_blast_test_out_summary"],self.id + ".csv"))
        self.genBlastReport("premrna_blast_control_out", join(dConfig["resources"]["premrna_blast_control_out_summary"],self.id + ".csv"))

    #@migrate 
    @updateResources
    def runMrnaBlast(self) -> None:

        if not "mrna" in self.fileDirectories:
            self.getResource("mrna")

        
        experiment = mrnaBlast.MrnaBlastExperiment(self.fileDirectories["mrna"])
        experiment.runExperiment()
        self.genBlastReport("mrna_blast_test_out", join(dConfig["resources"]["mrna_blast_test_out_summary"],self.id + ".csv"))
        self.genBlastReport("mrna_blast_control_out", join(dConfig["resources"]["mrna_blast_control_out_summary"],self.id + ".csv"))
