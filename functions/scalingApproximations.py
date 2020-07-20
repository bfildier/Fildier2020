"""Module scalingApproximations

Contains
- functions to compute the scaling approximations as in OGorman & Schneider (2009)
"""

#---- Modules -----#
import sys,os
from math import *
import numpy as np
from thermoFunctions import *
import xarray as xr

##-- O'Gorman&Schneider (2009) scaling --##

## Vertical integral on pressure coordinates
def verticalPressureIntegral(pres,values=None,dvdp=None,levdim=0):

    """Arguments: np.ndarray's or dask.array's with identical dimensions
    Returns: @f[ \int x \frac{dp}{g}@f] from bottom to top of atmosphere.
    It is defined negative for positive values.
    """
    
    cn = np
    
    dp = cn.diff(pres,axis=levdim) 

    val_prod = dp.copy()
    if values is not None:
        if values.__class__ in (np.ndarray,xr.core.dataarray.DataArray):
            nlev = pres.shape[levdim]
            val_mids = (np.take(values,range(nlev-1),axis=levdim) \
            + np.take(values,range(1,nlev),axis=levdim))/2
            val_prod = dp*val_mids
        if values.__class__ == list:
            nlev = pres.shape[levdim]
            val_prod = dp.copy()
            for i in range(len(values)):
                val_mids = (np.take(values[i],range(nlev-1),axis=levdim) \
                + np.take(values[i],range(1,nlev),axis=levdim))/2
                val_prod = val_prod*val_mids
    
    dvdp_prod = 1
    if dvdp is not None:
        if dvdp.__class__ in (np.ndarray,xr.core.dataarray.DataArray):
            dvdp_prod = dvdp
        elif dvdp.__class__ == list:
            for i in range(len(dvdp)):
                dvdp_prod = dvdp_prod*dvdp[i]

    return cn.nansum(val_prod*dvdp_prod/gg,axis=levdim)


## Vertical integral on pressure coordinates
def verticalPressureIntegralTroposphere(pres,temp,values=None,dvdp=None,levdim=0):

    """Arguments: np.ndarray's or dask.array's with identical dimensions
    Returns: @f[ \int x \frac{dp}{g}@f] from bottom to top of atmosphere.
    It is defined negative for positive values.
    
    DOES NOT WORK YET WITH dvdp
    """
    
    cn = np
    
    pres_c, temp_c, values_c = cropProfiles(pres,temp,values,levdim=levdim)
    
    dp = cn.diff(pres_c,axis=levdim) 

    val_prod = dp.copy()
    if values_c is not None:
        if values_c.__class__ == np.ndarray:
            nlev = pres_c.shape[levdim]
            val_mids = (np.take(values_c,range(nlev-1),axis=levdim) \
            + np.take(values_c,range(1,nlev),axis=levdim))/2
            val_prod = dp*val_mids
        if values_c.__class__ == list:
            nlev = pres_c.shape[levdim]
            val_prod = dp.copy()
            for i in range(len(values_c)):
                val_mids = (np.take(values_c[i],range(nlev-1),axis=levdim) \
                + np.take(values_c[i],range(1,nlev),axis=levdim))/2
                val_prod = val_prod*val_mids
    
    dvdp_prod = 1
    if dvdp is not None:
        if dvdp.__class__ == np.ndarray:
            dvdp_prod = dvdp
        elif dvdp.__class__ == list:
            for i in range(len(dvdp)):
                dvdp_prod = dvdp_prod*dvdp[i]

    return cn.nansum(val_prod*dvdp_prod/gg,axis=levdim)


## Find cold-point tropopause level
def tropopauseIndex(temp,levdim=0):
    
    cn = np

    return cn.argmin(temp,axis=levdim)
    # (or return when temperature first stops decreasing)
    # return temp1D.size - cn.argmax(cn.diff(cn.flipud(temp1D)) > 0)

## Find bottom index
def bottomIndex(pres,levdim=0):

    cn = np

    return cn.argmax(pres,axis=levdim)

## Crop profiles between surface and tropopause
def cropProfiles(pres,temp,values=None,levdim=0):

    cn = np

    # Replace 1D-data values v with nans outside slice s
    def cropVector(v,s):
        v_c = v.copy()
        v_c[:] = np.nan
        v_c[s] = v[s]
        return v_c
    # Apply cropVector function to all dimensions
    def cropArray(v,i_bot,i_trop,levdim):
        v_c = np.moveaxis(v,levdim,-1).copy()
        vshape = v_c.shape[:-1]
        for ind in cn.ndindex(vshape):
            i_min = min(i_bot[ind],i_trop[ind])
            i_max = max(i_bot[ind],i_trop[ind])
            s = slice(i_min,i_max+1)
            v_c[ind] = cropVector(v_c[ind],s)
        return np.moveaxis(v_c,-1,levdim)

    i_bot = bottomIndex(pres,levdim=levdim)
    i_trop = tropopauseIndex(temp,levdim=levdim)
    # Crop pressure and temperature arrays
    pres_c = cropArray(pres,i_bot,i_trop,levdim)
    temp_c = cropArray(temp,i_bot,i_trop,levdim)
    # Crop all other values arrays
    if values.__class__ == np.ndarray or values.__class__ == da.core.Array:
        values_c = cropArray(values,i_bot,i_trop,levdim)
#         print(values_c)
        return pres_c, temp_c, values_c
    elif values.__class__ == list:
        values_c = []
        for vals in values:
            values_c.append(cropArray(vals,i_bot,i_trop,levdim))
        return (pres_c,temp_c)+tuple(values_c)


## Compute O'Gorman & Schneider's scaling
def scalingOGS09_zcoord(wspeed,temp,pres,z,spechum=None,relhum=None,efficiency=1,temp_type='environment',parameter=1,
    levdim=0,ignore_qvstar=False):
    
    """Has not been tested for dask arrays."""

    cn = np

    if temp_type == 'environment':
        temp_profile = temp
#     elif temp_type == 'adiabat':
#         # Find out the vertical ordering of values
#         nlev = pres.shape[levdim]
#         i_b = bottomIndex(cn.moveaxis(pres,levdim,-1).ravel()[-nlev:])
#         reverse = False if i_b == 0 else True
#         # Compute T profile
#         temp_profile = moistAdiabatSimple(temp,pres,relhum=relhum,levdim=levdim,
#             reverse=reverse)
    elif temp_type == 'adiabat':
        # Find out the vertical ordering of values
        nlev = pres.shape[levdim]
        i_b = bottomIndex(cn.moveaxis(pres,levdim,-1).ravel()[-nlev:])
        reverse = False if i_b == 0 else True
        # Compute T profile
        temp_profile = moistAdiabatSimple(temp[i_b],pres,spechum=spechum,relhum=relhum,levdim=levdim,
            reverse=reverse)
#         print(temp_profile,temp_profile.shape)
    else:
        print("wrong type of temperature profile requested. Code this option.")
    # elif temp_type == 'parametric':
    #     p_ref = cn.take(pres,-1,axis=levdim)
    #     temp_profile = 

    pres_c, temp_c, wspeed_c, z_c = cropProfiles(pres,temp_profile,[wspeed,z],levdim=levdim)
    qvstar_c = saturationSpecificHumidity(temp_c,pres_c)

#     dqvstar_dp_c = cn.diff(qvstar_c,axis=levdim)/cn.diff(pres_c,axis=levdim)
    dqvstar_dz_c = cn.diff(qvstar_c,axis=levdim)/cn.diff(z_c,axis=levdim)
    
    if not ignore_qvstar:
        return -verticalPressureIntegral(pres_c,
                                         values=wspeed_c,
                                         dvdp=dqvstar_dz_c,
                                         levdim=levdim)
    else:
        return -verticalPressureIntegral(pres_c,
                                     values=wspeed_c,
                                     levdim=levdim)

def scalingOGS09_pcoord(omega,temp,pres,relhum=None,efficiency=1,temp_type='environment',parameter=1,
    levdim=0,ignore_qvstar=False):
    
    """Has not been tested for dask arrays."""

    cn = np

    if temp_type == 'environment':
        temp_profile = temp
    elif temp_type == 'adiabat':
        # Find out the vertical ordering of values
        nlev = pres.shape[levdim]
        i_b = bottomIndex(cn.moveaxis(pres,levdim,-1).ravel()[-nlev:])
        reverse = False if i_b == 0 else True
        # Compute T profile
        temp_profile = moistAdiabatSimple(temp,pres,relhum=relhum,levdim=levdim,
            reverse=reverse)
    else:
        print("wrong type of temperature profile requested. Code this option.")
    # elif temp_type == 'parametric':
    #     p_ref = cn.take(pres,-1,axis=levdim)
    #     temp_profile = 

    pres_c, temp_c, omega_c = cropProfiles(pres,temp_profile,omega,levdim=levdim)
    qvstar_c = saturationSpecificHumidity(temp_c,pres_c)
    dqvstar_dp_c = cn.diff(qvstar_c,axis=levdim)/cn.diff(pres_c,axis=levdim)
    
    if not ignore_qvstar:
        return -verticalPressureIntegral(pres_c,
                                         values=omega_c,
                                         dvdp=dqvstar_dp_c,
                                         levdim=levdim)
    else:
        return -verticalPressureIntegral(pres_c,
                                         values=omega_c,
                                         levdim=levdim)

## Compute O'Gorman & Schneider's scaling
def scalingAdditionL13_zcoord(wspeed,temp,pres,z,qv,relhum=None,efficiency=1,temp_type='environment',parameter=1,
    levdim=0,ignore_qvstar=False):
    
    """Has not been tested for dask arrays."""

    cn = np

    if temp_type == 'environment':
        temp_profile = temp
#     elif temp_type == 'adiabat':
#         # Find out the vertical ordering of values
#         nlev = pres.shape[levdim]
#         i_b = bottomIndex(cn.moveaxis(pres,levdim,-1).ravel()[-nlev:])
#         reverse = False if i_b == 0 else True
#         # Compute T profile
#         temp_profile = moistAdiabatSimple(temp,pres,relhum=relhum,levdim=levdim,
#             reverse=reverse)
    elif temp_type == 'adiabat':
        # Find out the vertical ordering of values
        nlev = pres.shape[levdim]
        i_b = bottomIndex(cn.moveaxis(pres,levdim,-1).ravel()[-nlev:])
        reverse = False if i_b == 0 else True
        # Compute T profile
        temp_profile = moistAdiabatSimple(temp[i_b],pres,spechum=qv,relhum=relhum,levdim=levdim,
            reverse=reverse)
    else:
        print("wrong type of temperature profile requested. Code this option.")
    # elif temp_type == 'parametric':
    #     p_ref = cn.take(pres,-1,axis=levdim)
    #     temp_profile = 

    pres_c, temp_c, wspeed_c, z_c, qv_c = cropProfiles(pres,temp_profile,[wspeed,z,qv],levdim=levdim)
    qvstar_c = saturationSpecificHumidity(temp_c,pres_c)

#     dqvstar_dp_c = cn.diff(qvstar_c,axis=levdim)/cn.diff(pres_c,axis=levdim)
#     dqvstar_dz_c = cn.diff(qvstar_c,axis=levdim)/cn.diff(z_c,axis=levdim)
    
    if not ignore_qvstar:
        return -verticalPressureIntegral(pres_c,
                                         values=[wspeed_c,qvstar_c-qv_c],
                                         levdim=levdim)
    else:
        return -verticalPressureIntegral(pres_c,
                                     values=wspeed_c,
                                     levdim=levdim)
    
def scalingL13_zcoord(wspeed,temp,pres,qv,temp_type='environment',
    levdim=0,ignore_qvstar=False):
    
    """Has not been tested for dask arrays."""

    cn = np

    if temp_type == 'environment':
        temp_profile = temp

    else:
        print("wrong type of temperature profile requested. Code this option.")

    pres_c, temp_c, wspeed_c, qv_c = cropProfiles(pres,temp_profile,[wspeed,qv],levdim=levdim)
    qvstar_c = saturationSpecificHumidity(temp_c,pres_c)
    prod_term = wspeed_c*(qvstar_c-qv_c)

    return -verticalPressureIntegral(pres_c,
                                     values=prod_term,
                                     levdim=levdim)
