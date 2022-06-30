"""
Developed for python3.8
justingreenblatt@github.com |last updated 09/05/2022
"""

from time import time
from os import environ
from os.path import join, isfile
from flask import Flask, render_template, url_for
from requests import packages

import myUtils
import logging
import sys

sys.path.insert(1, join(environ.get("HOME"), "blastWeb/bin"))
from ncbiData import ncbiData
from settings.logs import MAIN_LOG_PATH, MAIN_LOG_LEVEL
from settings.directories import NCBI_RESOURCE_LINKS_PATH
from ncbiScrape import getSpeciesLinks

if isfile(NCBI_RESOURCE_LINKS_PATH):
    h = open(NCBI_RESOURCE_LINKS_PATH)
    ncbiSpeciesLinks = json.load(h)
else:
    ncbiSpeciesLinks = getSpeciesLinks()

webGenomes = dict({a : list([ncbiData(b, c) for b,c in ncbiSpeciesLinks[a].items()]) for a in ncbiSpeciesLinks}) 
"""
app = Flask(__name__)

@app.route("/")
def template_test():
    return render_template("home.html", speciesList = workingGenomes)

@app.route("/<species>")
def speciesData(species):
    speciesObj = None
    genomeName = 'X'
    for i in workingGenomes:
        if i.name == species:
            speciesObj = i
            for r in rm:
                if r.species.lower().replace(' ','_') == i.id:
                    genomeName = r.genome
    data = speciesObj.histograms["rmSizeHist"].getDict()
    y1 = speciesObj.histograms2d["gwSizeHist2d"].getDict()["maxBinCount"]
    y2 = speciesObj.histograms2d["gwCSizeHist2d"].getDict()["maxBinCount"]
    hd = data["histogram"]
    bd = list([str(round(a)) for a in data["bins"]])
    return render_template("speciesData.html", s = speciesObj, histogram  = hd,
                           bins = bd, gName = genomeName, maxY = max(y1, y2))

if __name__ == "__main__":
    app.run(debug=True,host="0.0.0.0", port=80)
"""
