# RA, 2020-09-02


import logging

# # Note: The 'pickle' was made with this code:
# from pypath.core.network import Network
# op = Network.omnipath()
# op.save_to_pickle("op.dat")

# from pypath.core.network import Network
# op = Network.from_pickle("op.dat")



exit()

# FAILS

from pypath.share import settings
settings.setup(progressbars=True)


from pypath.core import intercell

ic = intercell.get_db()
print(ic)

exit()


# Issue

import os

from pypath.resources import data_formats
print(list(data_formats.ligand_receptor.keys()))

# WORKING
for (k, v) in data_formats.ligand_receptor.items():
    from pypath.legacy.main import PyPath
    pa = PyPath()
    print("Loading", k)
    pa.load_resources({k: v})
    print(pa.graph) #, file=open(os.devnull, mode='w'))


# OFFENSIVE
from pypath.legacy.main import PyPath
pa = PyPath()
for (k, v) in data_formats.ligand_receptor.items():
    if k not in ['lrdb', 'italk']:
        print("Loading", k)
        pa.load_resources({k: v})
        print(pa.graph, file=open(os.devnull, mode='w'))


exit()
