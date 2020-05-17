##
## Wrappers for pickling
##
import misopy
import os
import cPickle as pickle

def load_pickled_file(pickled_filename):
    if os.access(pickled_filename, os.F_OK):
        if os.name == 'nt':
            pickled_file = open(pickled_filename, 'rb') ### binary mode required for windows (EOFError otherwise at load)
        else:
            pickled_file = open(pickled_filename, 'r')
        loaded_obj = pickle.load(pickled_file)
        #loaded_obj = pickle.Unpickler(pickled_file)
        pickled_file.close()
        return loaded_obj
    return None

def write_pickled_file(obj_to_pickle, pickled_filename):
    if os.name == 'nt':
        pickled_file = file(pickled_filename, 'wb') ### binary mode required for windows (EOFError otherwise at load above) 
    else:
        pickled_file = file(pickled_filename, 'w')
    pickle.dump(obj_to_pickle, pickled_file, -1)
    #del pickled_file
    pickled_file.close()
