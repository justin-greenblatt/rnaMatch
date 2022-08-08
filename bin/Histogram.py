"""
Developed for python3.8
justinGreenblatt@github.com | last modified 09/05/2022
Histogram object for generating numpy histograms from data associated to an ensembleGenome object.
"""
import logging
from os.path import join
from numpy import histogram, histogram2d, array
from time import time
from settings.directories import HISTOGRAMS_FOLDER, HISTOGRAMS2D_FOLDER
import json
from sys import exit

#logging.basicConfig(filename= HISTOGRAM_LOG_PATH, level=HISTOGRAM_LOG_LEVEL)

class Histogram:
    """
    genome: Ensemble genome object
    dataKey: Key of for the genome object link dictionary that defines the source data for the histogram
    histkey: Key for the genome object histogram dictionary that will be used to store the relative path of the histogram generated here
    getDataFunc: lambda function that recieves a line of the data source file and extracts the desired dataPoint from it
    bins: number of bins used to generate histogram
    minMax: Array with to numeric values that define the minimum and maximum range of generated histogram
    skip = used to skip header of source data file. Numerical value equal to number of lines to skip
    """

    def __init__(self, genome, dataKey, histKey, getDataFunc, bins, minMax, skip=1):
        #logging.basicConfig(filename= HISTOGRAM_LOG_PATH, level=HISTOGRAM_LOG_LEVEL)
        #logging.info("Generating 1 dimensional histogram {} for {}|{}".format(histKey, genome.name, genome.species))

        #Searching for previously computed histogram for this genome
        if not histKey in genome.fileDirectories:
            #logging.debug("Generating 1D Histogram \t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(genome.name, dataKey, histKey, getDataFunc, bins, minMax, skip, time()))
            #load data source file
            dataHandler = open(genome.fileDirectories[dataKey])

            #skip headers
            for i in range(skip):
                dataHandler.readline()
            x = []
            for n,i in enumerate(dataHandler):
                #Get data points from source
                if n >= skip:
                    try:
                        x.append(getDataFunc(i))
                    except Exception as e:
                        print(e)
                        #logging.error("{}-{}\n{}".format(n,i,e))
                        exit(0)
            #Generate histogram with numpy.histogram function
            hist = histogram(x,bins,range=minMax)
            #Export data
            outPath = join(HISTOGRAMS_FOLDER, "{}_{}.json".format(genome.id, histKey))
            outHandler = open(outPath, 'w') 
            histDict = {"histogram":hist[0].tolist(), "bins": hist[1].tolist()}
            json.dump(histDict, outHandler)
            dataHandler.close()
            outHandler.close()
            #Set path in genome object to exported histogram data
            self.path =  outPath
            genome.fileDirectories[histKey] = outPath
            # In the case the data already exists just link genomoe object to path
        else:
            self.path = genome.fileDirectories[histKey]
            #logging.debug("Histogram exists at {} retrieving it \t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.path, genome.name, dataKey, histKey, getDataFunc, bins, minMax, skip, time()))


    def getDict(self):
        #loading histogram
        handler = open(self.path)
        jsonData = handler.read()
        handler.close()
        return json.loads(jsonData)



class Histogram2d(Histogram):
    """
    Generating 2 dimensional histogram
    genome: Ensemble genome object.
    dataKey: Key of for the genome object link dictionary that defines the source data for the histogram.
    histkey: Key for the genome object histogram dictionary that will be used to store the relative path of the histogram generated here.
    getDataFunc: lambda function that recieves a line of the data source file and extracts the desired ordered Pair (x,y) from it.
    xyBins: number of bins used to generate histogram in on the X dimension and Y dimension (xBins, yBins).
    xyMinMax: Array with to numeric values that define the minimum and maximum range of generated histogram ((xMin,xMax),(yMin, yMax)).
    skip = used to skip header of source data file. Numerical value equal to number of lines to skip.
    TIMESTAMP
    """
    def __init__(self, genome, dataKey, histKey, getDataFunc, xyBins, xyMinMax, skip=1):
        #logging.basicConfig(filename= HISTOGRAM2D_LOG_PATH, level=HISTOGRAM2D_LOG_LEVEL)
        #logging.info("Generating 2 dimensional histogram {} for {}|{}".format(histKey, genome.name, genome.species))
        #if histogram not generated yet
        if not histKey in genome.fileDirectories:
            #logging.debug("Generating 2D Histogram \t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(genome.name, dataKey, histKey, getDataFunc, xyBins, xyMinMax, skip, time()))
            dataHandler = open(genome.fileDirectories[dataKey])
            #filter data point (x,y) from data source lines
            xyPair = []
            for n,i in enumerate(dataHandler):
                if n >= skip:
                    try:
                        xyPair.append(getDataFunc(i))
                    except Exception as e:
                        print(e)
                        #logging.error("{}{}\n{}".format(n,i,e))
                        exit(0)
            dataHandler.close()
            #separate to two arrays the ordered pair
            x = list([p[0] for p in xyPair])
            y = list([p[1] for p in xyPair])
            hist = histogram2d(x, y, xyBins, xyMinMax)
            maxBinCount = hist[0].max()
            #export data
            outPath = join(HISTOGRAMS2D_FOLDER, "{}_{}.json".format(genome.id, histKey))
            outHandler = open(outPath,'w')
            histDict = {"histogram":hist[0].tolist(), "xbins": hist[1].tolist(), "ybins": hist[2].tolist(), "maxBinCount": maxBinCount}
            json.dump(histDict, outHandler)
            outHandler.close()
            dataHandler.close()
            genome.fileDirectories[histKey] = outPath
            self.path = outPath
        else:
            self.path =  genome.fileDirectories[histKey]
            #logging.debug("2D Histogram exists at {} retrieving it \t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.path, genome.name, dataKey, histKey, getDataFunc, xyBins, xyMinMax, skip, time()))
