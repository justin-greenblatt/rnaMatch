"""
Developed on python3.8
justingreenblatt@github.com | last updated 09/05/2022
The class has the functionality of retrieving data, holding paths to that data, generating reports, histograms and summaries of the data.

It contains all the information taken from genomes in the ensembl FTP site https://www.ensembl.org/info/data/ftp/index.html.
It contains over 305 genomes. A few of these are variants of the same species.
The constructer of the class recieves a species name that is on The ensembl FTP server. This information is previously scraped by
getEnsembleGenomes(). We call this species name a preLink cause it is an integral part of the FTP urls and has an all lower cap
and underline syntax. The class scrapes for all data links associated with the species in the ensenmbl FTP server and saves that in
a dictionary.
"""
import gzip
from subprocess import Popen, PIPE
from collections import Counter
from os import chdir, getcwd, remove, listdir
from os.path import join, isfile
from requests import get
from re import findall
from Bio import SeqIO
from settings.regularExpressions import ENSEMBLE_FTP_REGEX_GET_SPECIES
from settings.directories import GENOME_FOLDER, GENOME_WALK_FOLDER, GENOME_WALK_PATH, ENSEMBL_HTML_PATH, GTF_FOLDER, REPEAT_MASK_FOLDER, HISTOGRAMS_FOLDER, FILE_DIRECTORIES, LOG_FOLDER, HISTOGRAMS2D_FOLDER
from settings.links import ENSEMBL_DATA_LINK_PREFIX, ENSEMBL_FTP_LINK, ENSEMBL_LINK_TYPES
from settings.logs import ENSEMBL_GENOME_LOG_PATH, ENSEMBL_GENOME_LOG_LEVEL
from myUtils import downloadFromURL
from sys import exit
from numpy import histogram, histogram2d, arange
import json
import logging

logging.basicConfig(filename=ENSEMBL_GENOME_LOG_PATH, level = ENSEMBL_GENOME_LEVEL)

class ensemblGenome:

    def __init__(self, preLink, name, species):
        logging.info("Generating ensemblGenome object for {}|{}".format(name, species))
        logging.debug("Creating genome object: name species prelink [ {}, {}, {}]".format(name,species,preLink) 
        self.name = name
        self.preLink = preLink.lower()
        self.species = species
        self.id = species.lower().replace(" ","_")

        #These are all the datatypes that ensembl FTP offers for its species. This table holds information for inferring the url links to retrieve this data.
        dataTypes = ENSEMBL_LINK_TYPES
        self.linkDict = {a[0]: join(ENSEMBL_DATA_LINK_PREFIX, a[1], self.preLink, a[2]) for a in dataTypes}
        self.fileDirectories = {}
        logging.debug("ftp links related to genome object: dictionaryOfLinks [{}]".format(str(self.linkDict))

        #Listing directories in the program - This is inneficient. Should be done separrately
        genomeFiles = listdir(GENOME_FOLDER)
        gtfFiles = listdir(GTF_FOLDER)
        genomeWalkFiles = listdir(GENOME_WALK_FOLDER)
        repeatMaskFiles = listdir(REPEAT_MASK_FOLDER)
        histogramFiles = listdir(HISTOGRAMS_FOLDER)
        histogram2DFiles = listdir(HISTOGRAMS2D_FOLDER)

        for i in genomeFiles:
            if self.id in i.lower().replace(" ","_"):
                self.fileDirectories["dna"] = join(GENOME_FOLDER, i)
        for i in gtfFiles:
            if self.id in i.lower().replace(" ","_"):
                self.fileDirectories["gtfAnnotation"] = join(GTF_FOLDER, i)
        for i in genomeWalkFiles:
            if self.id in i.lower().replace(" ","_") and "control" not in i:
                self.fileDirectories["genomeWalk"] = join(GENOME_WALK_FOLDER, i)

        for i in genomeWalkFiles:
            if self.id in i.lower().replace(" ","_") and "control" in i:
                self.fileDirectories["genomeWalkControl"] = join(GENOME_WALK_FOLDER, i)

        for i in repeatMaskFiles:
            if self.id in i.lower().replace(" ","_"):
                self.fileDirectories["repeatMask"] = join(REPEAT_MASK_FOLDER, i)

        speciesHistograms = list([a for a in histogramFiles if self.id in a.lower().replace(" ","_")])
        for i in speciesHistograms:
            histName = i.split('/')[-1].split(".")[0].split("_")[-1]
            self.fileDirectories[histName] = join(HISTOGRAMS_FOLDER, i)

        speciesHistograms2D = list([a for a in histogram2DFiles if self.id in a.lower().replace(" ","_")])
        for i in speciesHistograms2D:
            histName = i.split('/')[-1].split(".")[0].split("_")[-1]
            self.fileDirectories[histName] = join(HISTOGRAMS2D_FOLDER, i)

        logging.debug("local data paths related to genome object: dictionaryOfPaths [{}]".format(str(self.fileDirectories))

    def getGenome(self):
    """Look for fasta genome associated to this object. If does not exist. Download genome from ensembl FTP 
    """

        if not "dna" in self.fileDirectories:

            getLink = self.linkDict["dna"]
            raw = get(getLink,verify=False).text
            genomeLink = join(getLink, findall("=\"(.*?sm.toplevel.fa.gz)",raw)[0])
            genomeAlreadyDownloaded = False
            old = getcwd()
            chdir(folder)
            self.fileDirectories["dna"] = join(GENOME_FOLDER,downloadFromURL(genomeLink))
            chdir(old)
            logging.debug("Genome Downloaded: name species directory [{},{},{}]".format(self.name, self.species, self.fileDirectories["dna"]))
        else:
            logging.debug("getGenome called but is Already Downloaded: name species directory [{},{},{}]".format(self.name, self.species, self.fileDirectories["dna"]))



    def generateHistograms(self):
    """Generate histograms of data associated to this object.
    """

         logging.info("generating histograms for {}|{}".format(self.name, self.species))

         self.histograms = {}
         self.histograms["gwSizeHist"] = Histogram(self, "genomeWalk", "gwSizeHist", lambda x: int(x.split(',')[-2]), 1000, (float(0),float(1000)))
         self.histograms["gwCSizeHist"] = Histogram(self,"genomeWalkControl", "gwCSizeHist", lambda x: int(x.split(',')[-2]), 1000, (0,1000))
         if "mSizeHist" in self.fileDirectories:
             self.fileDirectories["rmSizeHist"] = self.fileDirectories["mSizeHist"]
             self.fileDirectories.pop("mSizeHist")
         self.histograms["rmSizeHist"] = Histogram(self, "repeatMask", "rmSizeHist", lambda x: (lambda z: int(z[1]) - int(z[0]))([a for a in x.split() if a][5:7]), 1000, (0, 1000), 3)

         self.histograms2D = {}
         self.histograms2D["gwSizeHist2d"] = Histogram2d(self, "genomeWalk", "gwSizeHist2d", lambda x: int(x.split(',')[-2]), lambda x: float(x.split(',')[-1]), (1000,35), ((0, 1000), (0.65, 1)),skip = 3)
         self.histograms2D["gwCSizeHist2d"] = Histogram2d(self, "genomeWalk", "gwCSizeHist2d", lambda x: int(x.split(',')[-2]), lambda x: float(x.split(',')[-1]), (1000,35), ((0, 1000), (0.65, 1)), skip = 3)

         logging.debug("generated histograsm for ensemble object: name species histograms histograms2d [{},{},{},{}]".format(self.name, self.species, str(self.histograms), str(self.histograms2d)))

    def getGtf(self):

        logging.info("Getting gtf annontations for {}|{}".format(self.name, self.species))
        if not "gtfAnnotation" in self.fileDirectories:
            getLink = self.linkDict["gtfAnnotation"]
            raw = get(getLink,verify=False).text
            gtfLink = None
            try:
                gtfLink = findall("=\"(.*?\\.chr\\.gtf\\.gz)",raw)[0]
            except:
                gtfLink = findall("=\"(.*?\\.gtf\\.gz)",raw)[0]
            old = getcwd()
            chdir(GTF_FOLDER)
            self.fileDirectories["gtfAnnotation"] = join(GTF_FOLDER, downloadFromURL(join(getLink, gtfLink)))
            chdir(old)
            logging.debug("gtfAnnotations Downloaded: name species directory [{},{},{}]".format(self.name, self.species, self.fileDirectories["gtfAnnotation"]))
        else:
            logging.debug("getGtf called but is Already Downloaded: name species directory [{},{},{}]".format(self.name, self.species, self.fileDirectories["gtfAnnotation"]))


    def runGenomeWalk(self):
        """Run genomeWalk.py for this object. Download input gtf and genome files if they are not downloaded yet. Then create process
        """
        logging.info("Running genomeWalk for {}|{}".format(self.name, self.species))

        blastOutFile = join(GENOME_WALK_FOLDER, self.preLink + ".gw")
        blastControlOutFile = join(GENOME_WALK_FOLDER, self.preLink + "_control.gw")

        if not blastOutFile in self.fileDirectories:
            if not "gtfAnnotation" in self.fileDirectories:
                try:
                    self.getGtf()
                except:
                    logging.error("getGtf failed for {}|{}".format(self.name, self.species))

            if not "dna" in self.fileDirectories:
                try:
                    self.getGenome()
                except:
                    logging.error("getGenome failed for {}|{}".format(self.name, self.species))

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
                 logging.error("failed to run genomeWalk for {}|{}".format(self.name, self.species))

            genomeWalkFiles = listdir(GENOME_WALK_FOLDER)

            for i in genomeWalkFiles:
                if self.id in i.lower().replace(" ","_") and "control" not in i:
                    self.fileDirectories["genomeWalk"] = join(GENOME_WALK_FOLDER, i)

            for i in genomeWalkFiles:
                if self.id in i.lower().replace(" ","_") and "control" in i:
                    self.fileDirectories["genomeWalkControl"] = join(GENOME_WALK_FOLDER, i)
            logging.debug("Added paths to fileDirectories: genomeWalkPath genomeWalkControlPath [{},{}]".format(self.fileDirectories["genomeWalk"], self.fileDirectories["genomeWalkControl"]))

def getEnsemblGenomes():
"""This function generates ensemblGenome objects  for all species on the ensembl ftp page ensembl.org/info/data/ftp/index.html
"""
    ensembleHtmlData = get(ENSEMBL_FTP_LINK, verify=False).text
    return list([ensemblGenome(*a)
            for a in findall(ENSEMBLE_FTP_REGEX_GET_SPECIES, ensembleHtmlData)])
