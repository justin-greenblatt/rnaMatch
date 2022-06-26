from os.path import join
import re, requests, json
import logging
from time import time


def getSpeciesLinks():
    #logging.basicConfig(filename=NCBI_LOG_PATH, level = NCBI_LOG_LEVEL)
    #dissable warnings
    requests.packages.urllib3.disable_warnings(requests.packages.urllib3.exceptions.InsecureRequestWarning)
    def getLinks(link):

        return re.findall(r'href=\"(.*?)\"', requests.get(link, verify = False).text)


    GFF_COLUMNS = ['seqId', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'atributes']
    GFF_SKIPROWS = 8

    NCBI_LINK_START = "https://ftp.ncbi.nih.gov/genomes/refseq"
    NCBI_LINK_END = "annotation_releases/current"
    NCBI_EUKARYOTS = [
                      "vertebrate_mammalian",
                      "vertebrate_other",
                      "invertebrate",
                      "plant",
                      "protozoa",
                      "fungi"]
    #Find All eukaryot species in ncbi refseq
    links = dict({e :list([ a.strip('/') for a in getLinks(join(NCBI_LINK_START, e)) if a.endswith('/') and not a.startswith('/')]) for e in NCBI_EUKARYOTS})
    #Annotation Release Links
    arl = {a: {b: getLinks(join(NCBI_LINK_START,a,b,NCBI_LINK_END)) for b in links[a]} for a in links}
    #Filter Relevant Links
    filtered = {a: {b:[ c for c in arl[a][b] if len(c) > 5 and not c.startswith("mailto") and c.startswith("GCF")]for b in arl[a]} for a in arl}
    #drop species with no annotations
    dropEmpty = {a: {b: filtered[a][b] for b in filtered[a] if len(filtered[a][b]) > 0} for a in filtered}
    newLinks = {a: {b + '/' + dropEmpty[a][b][0] : {c.split('/')[-1] : join(NCBI_LINK_START,a,b,NCBI_LINK_END,dropEmpty[a][b][0],c) for c in getLinks(join(NCBI_LINK_START,a,b,NCBI_LINK_END,dropEmpty[a][b][0]))}  for b in dropEmpty[a]} for a in dropEmpty}


    return newLinks

if __name__ == "__main__":
    print("Downloading Links")
    h = open("ncbiLinks4.json",'w')
    l = dict(getSpeciesLinks())
    json.dump(l,h, indent = 2)
