from time import time
import sys
from os import environ
from os.path import join
from flask import Flask, render_template, url_for
from settings.directories import LOG_FILE

sys.path.insert(1, join(environ.get("HOME"), "blastWeb/bin"))

import ensemblData
import myUtils
import repeatMaskData
import logging

logging.basicConfig(filename=LOG_FILE, filemode='w', level=logging.INFO)

logging.info("Scraping data links from ensemble rest API and generating ensemble genome Objects")
genomes = ensemblData.getEnsemblGenomes()
logging.info("Scraping data links from repeatMask.org of program output")
rm = repeatMaskData.getAllData()
logging.info("Gettint intersection of species with data in ensembl and repeatMask.org")
webGenomes = list([g for g in genomes if g.species.lower() in [n.species.lower() for n in rm] and not "-" in g.name])
for g in webGenomes:
    try:
        g.generateHistograms()
    except:
        logging.error("Failed to generate histograms for {} | {}".format(g.species, g.name))
workingGenomes = list([a for a in webGenomes if a.histograms])

app = Flask(__name__)

@app.route("/")
def template_test():
    return render_template("home.html", speciesList = workingGenomes)


@app.route("/<species>")
def speciesData(species):
    speciesObj = None
    for i in workingGenomes:
        if i.name == species:
            speciesObj = i
    data = speciesObj.histograms["rmSizeHist"].getDict()
    hd = data["histogram"]
    bd = list([str(round(a)) for a in data["bins"]])
    return render_template("speciesData.html", s = speciesObj, histogram  = hd, bins = bd)

if __name__ == "__main__":
    app.run(debug=True,host="0.0.0.0", port=80)
