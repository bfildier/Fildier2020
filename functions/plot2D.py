
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt

## Plot vertical data (transect, or vertical profiles over time)
def subplotVerticalData(ax,x,y,Z,cmap=plt.cm.seismic,vmin=None,vmax=None,cbar=True):
    
    """Arguments:
        - x and y are coordinate values
        - Z are the data values flattened"""
    
    xs0,ys0 = np.meshgrid(x,y)
    xmin = min(x[0],x[-1])
    xmax = max(x[0],x[-1])
    ymin = min(y[0],y[-1])
    ymax = max(y[0],y[-1])
    X0 = np.vstack([xs0.flatten(),ys0.flatten()]).T
  
    extent = (xmax,xmin,ymin,ymax)

    # New coordinates
    xs,ys = np.meshgrid(np.linspace(xmin,xmax,num=len(x)),
                        np.linspace(ymin,ymax,num=len(y)))
    X = np.vstack([xs.flatten(),ys.flatten()]).T

    # New values
    resampled = griddata(X0,Z,X, method='cubic')
    resampled_2D = np.reshape(resampled,xs.shape)

    im = ax.imshow(np.flipud(resampled_2D),extent=extent,interpolation='bilinear',
                   cmap=cmap,vmin=vmin,vmax=vmax,aspect='auto',origin='upper')
    if cbar:
        plt.colorbar(im,ax=ax)
    
    if y[0] > y[-1]:
        plt.gca().invert_yaxis()