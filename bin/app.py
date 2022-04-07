from flask import Flask, render_template, url_for
import settings.directories

class sTest:
    def __init__(self,name):
        self.name = name
        self.species = name*3

app = Flask(__name__)

data = list([sTest(str(a)) for a in range(10)])

@app.route("/")
def template_test():
    return render_template("home.html", speciesList = data)

@app.route("/<species>")
def speciesData(species):
    speciesObj = None
    for i in data:
        if i.name == species:
            speciesObj = i

    return render_template("speciesData.html", s = speciesObj)
    
if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=80)
