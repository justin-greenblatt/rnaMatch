from flask import render_template, url_for
from ExonTrack import ExonTrack
class CircosPlot:

    def __init__(self, name):

        self.name = name
        self.trackList = []
        self.trackList.append(ExonTrack("ARC03", 150, 200))
        self.trackList.append(ExonTrack("ARC04", 250, 270))
        self.divName = self.name + "Div"
        self.objName = self.name + "Obj"


    def renderPlot(self):
        return render_template("test2.html",circosConfig = self)
