# RA, 2020-06-12

"""
Download the M. jannaschii genome.
The first archaeon to have its complete genome sequenced.

Refs:
    https://en.wikipedia.org/wiki/Methanocaldococcus_jannaschii
    https://science.sciencemag.org/content/273/5278/1058.long
    https://www.genome.jp/kegg-bin/show_organism?org=mja
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1
    https://www.ncbi.nlm.nih.gov/nuccore/L77117.1
"""

import wget
import hashlib
import pandas as pd

from pathlib import Path
from datetime import datetime, timezone

PARAM = {
    'fasta': "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1/GCA_000091665.1_ASM9166v1_genomic.fna.gz",
    'readme': "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1/README.txt",
    'md5': "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1/md5checksums.txt",

    'local': Path(__file__).parent.parent / "20200608-Downloads/genomes/M.jannaschii/UV/",
}


def download(url: str):
    out = PARAM['local'] / Path(url).name
    out.parent.mkdir(exist_ok=True, parents=True)

    if out.exists():
        return out

    wget.download(str(url).format(), str(out))

    with open(str(out) + "_meta.txt", 'w') as fd:
        print("Downloaded on {} from".format(datetime.now(tz=timezone.utc)), file=fd)
        print(url, file=fd)

    return out


def download_all():
    return {
        k: download(PARAM[k])
        for k in ['fasta', 'readme', 'md5']
    }


def check_md5(files):
    md5 = pd.read_csv(files['md5'], header=None, sep=' ').dropna(axis=1).set_index(0).squeeze()
    md5 = md5.apply(lambda f: Path(f).name)
    md5 = pd.Series(index=md5, data=md5.index)

    for (k, f) in files.items():
        if (f.name in md5):
            match = (md5[f.name] == hashlib.md5(f.open('rb').read()).hexdigest())
            match = ("OK" if match else "Failed")
            print(F"md5 check {f.name}: {match}")
        else:
            print(F"No md5 check for {f.name}")

    return files


def main():
    files = check_md5(download_all())
    print(files)

if __name__ == '__main__':
    main()
