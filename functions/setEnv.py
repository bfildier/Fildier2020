"""Module setEnv

Contains functions to define environment directories and other environment variables
"""

#---- modules ----#

import socket, os

hostname = socket.gethostname()
model = 'SAM6.10.10_EDMF'

## Test whether we are running this code on a laptop or a cluster (NERSC)
def isLaptop():

    return hostname in ['jollyjumper','tornado','clarity'] or "ucbvpn" in hostname

## Return mathine name
def getCurrentMachine():
    
    if isLaptop():
        return hostname
    else:
        return "coriknl"

## Define dataroot
def getDataroot():

    if isLaptop():
        ## Parent directory to load data
        return "/Users/bfildier/Data/simulations"
    else:
        return "/global/cscratch1/sd/bfildier"

## Define archivedir
def getArchivedir(model=model,machine='current'):
    
    dataroot = getDataroot()
    if machine == 'current':
        machine = getCurrentMachine()
    archivedir = os.path.join(dataroot,model,'archive',machine)

    return archivedir

## Define model directory
def getModelDir(model=model,machine='current'):

    if isLaptop():
        modeldir=os.path.join('/Users/bfildier/Code/',model)
    else:
        modeldir=os.path.join('/global/homes/b/bfildier/code',model)

    return modeldir

