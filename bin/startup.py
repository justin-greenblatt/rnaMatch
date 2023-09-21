from subprocess import Popen, PIPE
import os, json
from settings import dConfig as homelessConfig
from settings import configsPath
from ncbiScrape import getSpeciesLinks
#configure home Directorie
homelessConfig["common"]["HOME_DIR"] = os.environ.get("HOME")
dConfigOut = open('/' + str(os.path.join(*os.path.realpath(__file__).split('/')[:-1], "settings", "directories.ini")), 'w')
homelessConfig.write(dConfigOut)
dConfigOut.close()

from settings import sConfig, dConfig, pConfig

def runCommand(command):
    p = Popen(command, stdout = PIPE, stdin = PIPE)
    print(" ".join(command))
    p.wait()
    print("finished")
    return p.communicate()

for p in sConfig["apt"].values():
    runCommand(["sudo", "apt", "install", p, "-y"])
for q in sConfig["pip"].values():
    runCommand(["sudo", "pip3", "install", q, "--no-input"])

def createDir(dirName):
    if not os.path.isdir(dirName):
        os.mkdir(dirName)       
        runCommand(["chmod", "-R", "777", dirName])

createDir(dConfig["common"]["nfs_path"])
createDir(dConfig["common"]["data_path"])
#mountDiskCommand = ["sudo", "mount", sConfig["nfs"]["IP"], dConfig["common"]["nfs_path"]]
#c = runCommand(mountDiskCommand)
#print(c)
for f in dConfig["resourceFolders"].values():
    createDir(f)

for d in dConfig["resources"].values():
    createDir(d)

for n in dConfig["nfs"].values():
    createDir(n)

links = 0
for ln in dConfig["ncbiLinks"].values():
    createDir(ln)
    links += len(os.listdir(ln))

if links == 0:
    print("Downoading Links from ncbiFtp")
    scrapeLinks = getSpeciesLinks()
    for k in scrapeLinks:
        groupFolder = dConfig["ncbiLinks"][k]
        for s in scrapeLinks[k]:
            h = open(os.path.join(groupFolder,s.rstrip('/').replace('/','-') + '.json'), 'w')
            json.dump(scrapeLinks[k][s],h)
            h.close() 
