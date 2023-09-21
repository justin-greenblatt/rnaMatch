"""
Justin Greenblatt 27/06/2022
Function for retrieving all Eukaryot current annotation assemblies and resources associated to Them.
"""


from os.path import join
import re, requests, json
from time import strftime ,time
import logging

from settings.links import NCBI_LINK_START, NCBI_LINK_END, NCBI_EUKARYOTS

def getSpeciesLinks():

    #setting up logging
    #logging.basicConfig(filename=NCBI_SCRAPE_LOG_PATH, level = NCBI_SCRAPE_LOG_LEVEL)
    #dissable warnings
    requests.packages.urllib3.disable_warnings(requests.packages.urllib3.exceptions.InsecureRequestWarning)

    #logging.info("Starting webScrape of ncbi ftp. time = {}".format(TIME()))

    def getLinks(link):
        #Util function for extracting all links from an hml response. Used extensively in the list comprehension to follow.
        t0 = time()
        print("Get request to {} ".format(link))
        linksFound = re.findall(r'href=\"(.*?)\"', requests.get(link, verify = False).text)
        print("request took {} seconds".format(round(time() - t0),1))
        return linksFound

    #Find All eukaryot species in ncbi refseq
    links = dict({e :list([ a.strip('/') for a in getLinks(join(NCBI_LINK_START, e)) if a.endswith('/') and not a.startswith('/')]) for e in NCBI_EUKARYOTS})
    print("Found assemblies for the following groups:{}".format(str(dict({a:len(b) for a,b in links.items()}))))
 
    #Annotation Release Links
    arl = {a: {b: getLinks(join(NCBI_LINK_START,a,b,NCBI_LINK_END)) for b in links[a]} for a in links}

    #Filter Relevant Links
    filtered = {a: {b:[ c for c in arl[a][b] if len(c) > 5 and not c.startswith("mailto") and c.startswith("GCF")]for b in arl[a]} for a in arl}

    #drop species with no annotations
    dropEmpty = {a: {b: filtered[a][b] for b in filtered[a] if len(filtered[a][b]) > 0} for a in filtered}
    newLinks = {a: {b + '/' + dropEmpty[a][b][0] : {c.split('/')[-1] : join(NCBI_LINK_START,a,b,NCBI_LINK_END,dropEmpty[a][b][0],c) for c in getLinks(join(NCBI_LINK_START,a,b,NCBI_LINK_END,dropEmpty[a][b][0]))}  for b in dropEmpty[a]} for a in dropEmpty}


    return newLinks

#Test for this unit

if __name__ == "__main__":
    print("Downloading Links")
    h = open("ncbiLinks.json",'w')
    l = dict(getSpeciesLinks())
    json.dump(l,h, indent = 2)

