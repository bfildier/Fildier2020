"""Module dataFormat

Define parameters related to the use of dask arrays that are common across modules
"""

#---- Modules ----#
import numpy as np
import re
from math import modf

# ## Default compute action for dask arrays
# da_compute_default = False
# ## Size of dask array chunks
# chunks = 100000

#---- Functions ----#

## Get module corresponding to array type
def getArrayType(array):

	if array.__class__.__bases__[0] is np.ndarray or array.__class__ is np.ndarray:
		cn = np # stands for 'class name'
	# elif array.__class__ == da.core.Array:
	# 	cn = da
	else:
		print("[Error in daskOptions.getArrayType] Unvalid data type:", type(array))
		return 

	return cn

## Turn string representing a decimal number into float, as for CS in experiment name
def str2float(s):
    
    pattern_dec = re.compile("^[0-9]+d[0-9]+$")
    pattern_0 = re.compile("^0")
    if pattern_dec.match(s):
        return float(s.replace('d','.'))
    elif pattern_0.match(s):
        return float(s[0]+'.'+s[1:])
    else:
        return float(s)

## Reverse transformation
def float2str(f,ndec=4,sep=''):
    
    """f: float to convert
    ndec: number of digits"""
    
    if (0<f<1 or -1<f<0) and sep == '':
        stemp = "%%0.%df"%ndec
        s = (stemp%f).rstrip('0')
        return s.replace('.','')
    elif modf(f)[0] == 0:
        return str(int(f))
    else:
        f_dec,f_int = modf(f)
        s_int = float2str(f_int,ndec=ndec)
        s_dec = float2str(abs(f_dec),ndec=ndec)
        return s_int+sep+s_dec[1:]

## Extract information from experiment name
def getDayOrigin(expname):
    if re.match(".*\-b[0-9]+\-.*",expname):
        chunks = expname.split('-')
        bchunk = [c for c in chunks if c[0] == 'b'][0]
        return int(bchunk[1:])
    return 0