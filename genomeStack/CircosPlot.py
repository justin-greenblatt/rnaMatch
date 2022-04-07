from flask import render_template, url_for

class CircosPlot:

    def __init__(self, name, genomicRegion, templateDirectory = "./templates", staticDirectory = "./static", circosTemplate = "circosTemplate.html"):

        self.circosTemplate = circosTemplate
        self.templateDiretcory = templateDirectory
        self.static = staticDirectory
        self.genomicRegion = genomicRegion
        self.name = name
        self.genomeName = self.name + "Genome"
        self.chrom = self.genomicRegion.chrom
        self.start = self.genomicRegion.start
        self.end = self.genomicRegion.end
        self.chromSpacing = 0.04
        self.width = 900
        self.height = 600
        self.innerRadius = 246
        self.outerRadius = 250
        self.genomeFillColorList = ["#FF5733"]
        self.chromSize = self.end - self.start
        self.trackList = []
        self.trackList.append(ExonTrack("ARC01", 150, 200))
        self.trackList.append(ExonTrack("ARC02", 250, 270))
        self.divName = self.name + "Div"
        self.objName = self.name + "Obj"


    def renderPlot(self):
        return render_template("" + self.circosTemplate,circosConfig = self)

    def colors(self):
        return str(self.genomeFillColorList)
