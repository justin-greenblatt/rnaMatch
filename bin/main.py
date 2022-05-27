"""
Developed for python3.8
justingreenblatt@github.com |last updated 09/05/2022
"""
from time import time
import sys
from os import environ
from os.path import join
from flask import Flask, render_template, url_for
from settings.logs import MAIN_LOG_PATH, MAIN_LOG_LEVEL
from requests import packages
sys.path.insert(1, join(environ.get("HOME"), "blastWeb/bin"))
import ensemblData
import myUtils
import repeatMaskData
import logging

logging.basicConfig(filename = MAIN_LOG_PATH, level = MAIN_LOG_LEVEL)
packages.urllib3.disable_warnings(packages.urllib3.exceptions.InsecureRequestWarning)
logging.info("Scraping data links from ensemble rest API and generating ensemble genome Objects")
genomes = ensemblData.getEnsemblGenomes()
logging.info("Scraping data links from repeatMask.org of program output")
rm = repeatMaskData.getLinks()
logging.info("Gettint intersection of species with data in ensembl and repeatMask.org")
print(len(rm), "REPEAT MASK")
print(len(genomes), "ENSEMBL GENOMES")
webGenomes = list([g for g in genomes if g.species.lower() in [n.species.lower() for n in rm] and not "-" in g.name])
print(len(webGenomes), "WEBgENOMES")
for g in webGenomes:
        g.generateHistograms()

workingGenomes = list([a for a in webGenomes if a.histograms2d])
print(len(workingGenomes), "WORKINGGENOMES")
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
