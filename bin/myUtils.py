from requests import get
from subprocess import Popen
from os import listdir, getcwd
from shutil import unpack_archive
from os import listdir

def downloadFromURL(url, filename = False, decompress = False):
    if not filename:
        filename = url.split('/')[-1]
    h = open(filename, 'wb+')
    r = get(url, stream=True, allow_redirects=True, verify=False)
    data = r.content
    h = open(filename, 'wb+')
    h.write(data)
    h.close()

    if decompress:
        gunzip = Popen(["gunzip", filename])
        gunzip.wait()
        newName = ""

        for f in listdir('.'):
            if filename.startswith(f) and len(f) > len(newName):
                newName = f

    return filename
