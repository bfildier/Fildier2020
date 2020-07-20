


import pickle
import os


### Functions ###
def save2object(obj,filepath,overwrite=False):
    # Load preexisting results
    if os.path.exists(filepath) and obj.__class__ == dict and not overwrite:
        obj_old = pickle.load(open(filepath,'rb'))
        obj_old.update(obj)
        pickle.dump(obj_old,open(filepath,"wb"))
    else:
        pickle.dump(obj,open(filepath,"wb"))

