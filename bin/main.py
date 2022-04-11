from time import time
import sys
from os import environ
from os.path import join
from flask import Flask, render_template, url_for
from settings.directories import LOG_FOLDER

sys.path.insert(1, join(environ.get("HOME"), "GenomeBlastWalk/bin"))

import ensemblData
import myUtils
import repeatMaskData
import logging

logging.basicConfig(filename=LOG_FOLDER, filemode='w', level=logging.INFO)

logging.info("Scraping data links from ensemble rest API and generating ensemble genome Objects")
genomes = ensemblData.getEnsemblGenomes()
logging.info("Scraping data links from repeatMask.org of program output")
rm = repeatMaskData.getAllData()
logging.info("Gettint intersection of species with data in ensembl and repeatMask.org")
webGenomes = list([g for g in genomes if "repeatMask" in g.fileDirectories.keys() and "genomeWalk" in g.fileDirectories.keys() and not g.name.endswith(" ") and not g.name.endswith("-")])

for g in webGenomes:
    g.generateHistograms()

app = Flask(__name__)

@app.route("/")
def template_test():
    return render_template("home.html", speciesList = webGenomes)

@app.route("/<species>")
def speciesData(species):
    speciesObj = None
    for i in webGenomes:
        if i.name == species:
            speciesObj = i

    return render_template("speciesData.html", s = speciesObj)

#if __name__ == "__main__":
#    app.run(debug=True, host="0.0.0.0", port=80)
