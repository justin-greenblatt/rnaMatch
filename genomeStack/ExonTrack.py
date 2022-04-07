class ExonTrack:
    def __init__(self, name, innerRadius, outerRadius):
        class exon:
            def __init__(self, chrom, start, end, color, des):
                self.chrom = chrom
                self.start = start
                self.end = end
                self.color = color
                self.des = des
                self.js = "{chr: \"" + chrom + "\", start: \"" + start + "\", end: \"" + end + "\", color: \"" + color + "\", des: \"" + des + "\"}, "

        self.name = name
        self.innerRadius = innerRadius #smaller
        self.outerRadius = outerRadius #larger
        self.arcList = []
        for i in range(25):
            self.arcList.append(exon("I", str(50 + i*100), str(100 + i*100),"#CD8500",str(i)))