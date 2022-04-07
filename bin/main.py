from time import time
import sys
from os import environ
from os.path import join
from flask import Flask, render_template, url_for
import settings.directories

sys.path.insert(1, join(environ.get("HOME"), "GenomeBlastWalk/bin"))

import ensemblData
import myUtils
import repeatMaskData

genomes = ensemblData.getEnsemblGenomes()
rm = repeatMaskData.getAllData()
workingRm = list([g for g in [a for a in rm if a.species.lower() in [s.species.lower() for s in genomes]]])
workingGenomes = list([g for g in [a for a in genomes if a.species.lower() in [s.species.lower() for s in rm]]])
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
