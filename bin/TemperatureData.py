import pandas as pd
from settings.directories import ANAGE_DATA_FILE
from sklearn.cluster import KMeans
import numpy as np
from settings.colors import PASTEL5 as SHAPES_PALLETE
from settings.colors import COLORS5 as CLUSTERS_PALLETE 
from random import randint

class TemperatureData:

    #From a list of species, generate clusters through KM based on their temperature data.
    def __init__(self,speciesList, nClusters):
        #An object representing clusters to be past to home.html
        class Cluster:
            def __init__(self, dataFrame, name, color, opacity = 1, size = 10, idNumber = randint(0,100000), text = True):
                self.data = dataFrame
                self.name = name
                self.size = size
                self.id = idNumber
                if text:
                    self.textArray = self.data["sciName"].to_list()
                else:
                    self.textArray = list([' ' for a in self.data["sciName"]])

                self.xArray = list(round(a,2) for a in self.data["celcius"])
                self.yArray = list(round(a,2) for a in self.data["nOrder"])
        #An object to draw on the graph background. passed to home.html
        class Shape:
            def __init__(self, yStart, yEnd, color, name = False, opacity = 0.2):
                self.yStart = yStart
                self.yEnd = yEnd
                self.color = color
                self.opacity = opacity
        self.shapes = []
        self.clusters = []
        data = pd.read_table(ANAGE_DATA_FILE)

        #Generate ids from data entries
        data["speciesID"] = data["Genus"].apply(lambda x: x.lower()) + '_' + data["Species"]
        data["sciName"] = data["Genus"]  + ' ' + data["Species"]
        data["celcius"] = data["Temperature (K)"] - 273

        tempData = data[data["Temperature (K)"].notna()]
        classDict = {b:a for a,b in enumerate(list(tempData.Class.unique()))}
        inverseClassDict = {a:b for a,b in enumerate(list(tempData.Class.unique()))}
        classList = dict({a: [] for a,b in enumerate(list(tempData.Class.unique()))})
        tempData["nClass"] = tempData.Class.apply(lambda x: classDict[x])
        #cData = tempData.sort_values(by=['nClass'], inplace=True)
        
        orderDict = {b:a for a,b in enumerate(list(tempData.Order.unique()))}
        print(list(tempData.Order.unique()))
        self.yTicks = list(orderDict.keys())
        self.yTickValues = list(orderDict.values())
        tempData["nOrder"] = tempData.Order.apply(lambda x: orderDict[x])

        #outData = cData.sort_values(by=['nClass'], inplace=True)

        for a,b in zip(tempData.nClass, tempData.nOrder):
            classList[a].append(b)
        
        for m,n in classList.items():
            self.shapes.append(Shape(min(n), max(n), SHAPES_PALLETE[m],
                               name = inverseClassDict[m]))

        myData = tempData[tempData.speciesID.isin(speciesList)]
        km = KMeans(nClusters)
        km.fit(np.array(myData.celcius).reshape(-1, 1))
        myData["nCluster"] = km.labels_
        tempCluster = Cluster(tempData, "Temperature Data total:{}".format(tempData.shape[0]), "#777777", 0.5, idNumber = 0,  text = False)
        self.clusters.append(tempCluster) 
        myClusters = myData.nCluster.unique()
        for c in myClusters:
            clusterData = myData[myData.nCluster == c]
            
            newCluster =  Cluster(clusterData, "KM Cluster {} total:{}".format(c + 1, clusterData.shape[0]), CLUSTERS_PALLETE[c],idNumber = c + 1)
            self.clusters.append(newCluster)
