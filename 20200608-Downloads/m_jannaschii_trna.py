# RA, 2020-06-16

"""
Download the tRNA genomic sequences for M. jannaschii.
"""

import wget

from pathlib import Path
from datetime import datetime, timezone

ROOT_PATH = Path(__file__).parent

PARAM = {
    'urls': {
        # http://gtrnadb2009.ucsc.edu/Meth_jann/Meth_jann-align.html
        "http://gtrnadb2009.ucsc.edu/Meth_jann/methJann1-tRNAs.fa",
    },

    'local': ROOT_PATH / "genomes/M.jannaschii/trna/"
}


def download(url: str):
    assert url

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
    for url in PARAM['urls']:
        print(download(url))


if __name__ == '__main__':
    download_all()
