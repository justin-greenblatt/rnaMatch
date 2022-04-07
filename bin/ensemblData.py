import gzip
from subprocess import Popen, PIPE
from collections import Counter
from os import chdir, getcwd, remove, listdir
from os.path import join, isfile
from requests import get
from re import findall
from Bio import SeqIO
from settings.regularExpressions import ENSEMBLE_FTP_REGEX_GET_SPECIES
from settings.directories import GENOME_FOLDER, GENOME_WALK_FOLDER, GENOME_WALK_PATH, ENSEMBL_HTML_PATH, GTF_FOLDER, REPEAT_MASK_FOLDER, HISTOGRAMS_FOLDER
from settings.links import ENSEMBL_DATA_LINK_PREFIX, ENSEMBLE_FTP_LINK
from myUtils import downloadFromURL
from sys import exit
from numpy import histogram, histogram2d, arange
import json

class ensemblGenome:
    def __init__(self, preLink, name, species):
        self.name = name
        self.preLink = preLink.lower()
        self.species = species
        self.id = species.lower().replace(" ","_")
        dataTypes = (
                     ("dna", "fasta", "dna"),
                     ("cdna", "fasta", "cdna"),
                     ("ncrna", "fasta", "ncrna"),
                     ("protein", "fasta", "pep"),
                     ("gtfAnnotation", "gtf", ""),
                     ("gff3Annotation", "gff3", ""),
                     ("genebankAnnotatedSeq", "genebank", ""),
                     ("emblAnnotatedSeq", "embl", ""),
                     ("tsvAnnotation", "tsv", ""),
                     ("rdfAnnotation", "rdf", ""),
                     ("jsonAnnotation", "json", ""),
                     ("WholeDatabase", "mysql", ""),
                     ("gvfVariation", "gvf", ""),
                     ("vcfVariation", "vcf", ""),
                     ("gffRegulation", "regulation", ""),
                     ("regulationDataFiles", "data_files", ""),
                     ("bamBigWig", "bamcov","")
                     )

        self.linkDict = {a[0]: join(ENSEMBL_DATA_LINK_PREFIX, a[1], self.preLink, a[2]) for a in dataTypes}
        self.fileDirectories = {}

        genomeFiles = listdir(GENOME_FOLDER)
        gtfFiles = listdir(GTF_FOLDER)
        genomeWalkFiles = listdir(GENOME_WALK_FOLDER)
        repeatMaskFiles = listdir(REPEAT_MASK_FOLDER)
        histogramFiles = listdir(HISTOGRAMS_FOLDER)

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
            histName = i.lstrip(self.id + "_").rstrip(".json")
            self.fileDirectories[histName] = join(HISTOGRAMS_FOLDER, i)

    def getGenome(self, folder = GENOME_FOLDER):

        #Downloading ensemble Links for genome
        getLink = self.linkDict["dna"]
        raw = get(getLink,verify=False).text
        genomeLink = join(getLink, findall("=\"(.*?sm.toplevel.fa.gz)",raw)[0])

        genomeAlreadyDownloaded = False
        for f in listdir(GENOME_FOLDER):
            print("downloaded",join(GENOME_FOLDER, f))
            if self.preLink.lower() in f.lower():
                genomeAlreadyDownloaded = True
                self.fileDirectories["dna"] = join(GENOME_FOLDER,f)

        if isfile(join(GENOME_FOLDER, genomeLink.split('/')[-1])):
            genomeAlreadyDownloaded = True
            self.fileDirectories["dna"] = join(GENOME_FOLDER, genomeLink.split('/')[-1])

        if not genomeAlreadyDownloaded:
        #Download Genome
            old = getcwd()
            chdir(folder)
            self.fileDirectories["dna"] = join(GENOME_FOLDER,downloadFromURL(genomeLink))
            chdir(old)
            print("GenomeDownloaded")
        else:
            print("GENOME ALREADY DOWNLOADED")

    def generateHistograms(self):

         class Histogram:
              def __init__(self, genome, dataKey, histKey, getDataFunc, bins, minMax, skip=1):

                   if not histKey in genome.fileDirectories:
                        dataHandler = open(genome.fileDirectories[dataKey])
                        for i in range(skip):
                             dataHandler.readline()
                        x = list([getDataFunc(a) for a in dataHandler])
                        outDir = open(join(HISTOGRAMS_FOLDER, "{}_{}.json".format(genome.id, histKey)),'w')
                        print("GENERATING - ", genome.id, histKey)
                        hist = histogram(x,bins,range=minMax)

                        histDict = {"histogram":hist[0].tolist(), "bins": hist[1].tolist()}
                        json.dump(histDict, outDir)
                        dataHandler.close()
                        outDir.close()
                   else:
                        pass
                        print("loading - ", genome.fileDirectories[histKey])

         class Histogram2d:
              def __init__(self, genome, dataKey, histKey, getDataFuncX, getDataFuncY, XYbins, XYMinMax, skip=1):

                   if not histKey in genome.fileDirectories:

                        dataHandler = open(genome.fileDirectories[dataKey])
                        for i in range(skip):
                            dataHandler.readline()
                        x = list([getDataFuncX(a) for a in dataHandler])
                        dataHandler.close()

                        dataHandler = open(genome.fileDirectories[dataKey])
                        for i in range(skip):
                            dataHandler.readline()
                        y = list([getDataFuncY(a) for a in dataHandler])
                        dataHandler.close()

                        outDir = open(join(HISTOGRAMS_FOLDER, "{}_{}.json".format(genome.id, histKey)),'w')
                        hist = histogram(x, bins, (min, max))
                        histDict = {"histogram":hist[0].tolist(), "Xbins": hist[1].tolist(), "Ybins": hist[2].tolist()}
                        json.dumps(histDict, outDir)
                        outDir.close()

         self.histograms = {}
         self.histograms["gwSizeHist"] = Histogram(self, "genomeWalk", "gwSizeHist", lambda x: int(x.split(',')[-2]), 1000, (float(0),float(1000)))
         self.histograms["gwCSizeHist"] = Histogram(self,"genomeWalkControl", "gwCSizeHist", lambda x: int(x.split(',')[-2]), 1000, (0,1000))
         self.histograms["rmSizeHist"] = Histogram(self, "repeatMask", "rmSizeHist", lambda x: (lambda z: int(z[1]) - int(z[0]))([a for a in x.split() if a][5:7]), 1000, (0, 1000), 3)
         #self.histograms["gwSizeHist2d"] = Histogram2d(self, "genomeWalk", "gwSizeHist2d", lambda x: int(x.split(',')[-2]), lambda x: float(x.split(',')[-1]), (1000,35), ((0, 1000), (0.65, 1)))
         #self.histograms["gwCSizeHist2d"] = Histogram2d(self, "genomeWalk", "gwCSizeHist2d", lambda x: int(x.split(',')[-2]), lambda x: float(x.split(',')[-1]), (1000,35), ((0, 1000), (0.65, 1)))

    def parseGenome(self):

        self.chromossomes = {}
        print("Parsing",self.fileDirectories["dna"])
        genomeHandler = gzip.open(self.fileDirectories["dna"], "rt")
        genomeIterator = SeqIO.parse(genomeHandler, "fasta")
        self.baseCounter = Counter()
        self.genomeBaseData = Counter()
        self.softMaskSizes = Counter()
        for c in genomeIterator:
            print("Parsing", c.id)
            sequence = str(c.seq)
            self.chromossomes[c.id] = len(sequence)
            maskedSize = 0

            for n in sequence:
                self.genomeBaseData[n] +=1
                if n.islower():
                    maskedSize +=1
                else:
                    if maskedSize > 0:
                        self.softMaskSizes[maskedSize] +=1
                        maskedSize = 0

    def getGtf(self):
        getLink = self.linkDict["gtfAnnotation"]
        raw = get(getLink,verify=False).text
        gtfLink = None
        try:
            gtfLink = findall("=\"(.*?\\.chr\\.gtf\\.gz)",raw)[0]
        except:
            gtfLink = findall("=\"(.*?\\.gtf\\.gz)",raw)[0]

        genomeAlreadyDownloaded = False
        for f in listdir(GTF_FOLDER):
            print("downloaded",join(GTF_FOLDER, f))
            if self.preLink.lower() in f.lower():
                genomeAlreadyDownloaded = True
                self.fileDirectories["gtfAnnotation"] = join(GTF_FOLDER,f)

        if isfile(join(GTF_FOLDER, gtfLink.split('/')[-1])):
            genomeAlreadyDownloaded = True
            self.fileDirectories["gtfAnnotation"] = join(GTF_FOLDER, gtfLink.split('/')[-1])

        if not genomeAlreadyDownloaded:
            old = getcwd()
            chdir(GTF_FOLDER)
            self.fileDirectories["gtfAnnotation"] = join(GTF_FOLDER, downloadFromURL(join(getLink, gtfLink)))
            chdir(old)
        else:
            print("GTF FILE ALREADY DOWNLOADED")

    def runGenomeWalk(self):
        blastOutFile = join(GENOME_WALK_FOLDER, self.preLink + ".gw")
        blastControlOutFile = join(GENOME_WALK_FOLDER, self.preLink + "_control.gw")

        if not blastOutFile in self.fileDirectories:
            if not "gtfAnnotation" in self.fileDirectories:
                try:
                    self.getGtf()
                except:
                    print("error getting gtf for " + self.name)
                    exit(0)
            if not "dna" in self.fileDirectories:
                try:
                    self.getGenome()
                except:
                    print("error getting genome for " + self.name)
                    exit(0)
            try:
                 p = Popen(
                      ["python3", GENOME_WALK_PATH,
                      self.fileDirectories["dna"],
                      self.fileDirectories["gtfAnnotation"],
                      blastOutFile,
                      blastControlOutFile],
                      stdout = PIPE)
                 print("!!starting", self.species)
                 p.wait()
                 print("!!finished", self.species)
                 remove(self.fileDirectories["dna"])
                 remove(self.fileDirectories["gtfAnnotation"])

            except:
                pass
        self.fileDirectories.pop("gtfAnnotation")

def getEnsemblGenomes():
    h = open(ENSEMBL_HTML_PATH)
    ensembleHtmlData = h.read()
    h.close()
    return list([ensemblGenome(*a)
            for a in findall(ENSEMBLE_FTP_REGEX_GET_SPECIES, ensembleHtmlData)])
