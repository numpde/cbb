# RA, 2020-09-11

from os import devnull
from pypath.legacy.main import PyPath
from pypath.resources import data_formats

devnull = open(devnull, mode='w')

# print(list(data_formats.pathway.keys()))
# ['trip', 'spike', 'signalink3', 'guide2pharma', 'ca1', 'arn', 'nrf2', 'macrophage', 'death', 'pdz', 'signor', 'adhesome', 'icellnet', 'hpmr', 'cellphonedb', 'ramilowski2015', 'lrdb', 'baccin2019']

set1 = {"lrdb", "macrophage", "signalink3", "signor"}
set2 = set(data_formats.pathway.keys()) - set1

datasets = sorted(data_formats.pathway.items())

print("= Set 1 works =")
pa = PyPath()
for (k, v) in datasets:
    if k in set1:
        print(F"Loading {k}")
        pa.load_resources({k: v})
        print(pa.graph, file=devnull)

print("= Set 2 works =")
pa = PyPath()
for (k, v) in datasets:
    if k in set2:
        print(F"Loading {k}")
        pa.load_resources({k: v})
        print(pa.graph, file=devnull)

print("= Combination loads =")
pa = PyPath()
for (k, v) in datasets:
    print(F"Loading {k}")
    pa.load_resources({k: v})

print("= The graph fails =")
print(pa.graph)

