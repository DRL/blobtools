#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""usage: blobtools nodesdb             --nodes <NODES> --names <NAMES>
                                        [-h|--help]

    Options:
        -h --help                       show this
        --nodes <NODES>                 NCBI nodes.dmp file.
        --names <NAMES>                 NCBI names.dmp file.
"""

# from __future__ import division
from docopt import docopt

from os.path import basename, isfile, join, dirname, abspath
from sys import path
path.append(dirname(dirname(abspath(__file__))))

import lib.blobtools as blobtools
import lib.BtCore as BtCore
import lib.BtLog as BtLog
import lib.BtIO as BtIO

def main():
    args = docopt(__doc__)
    names_f = args['--names']
    nodes_f = args['--nodes']

    # Parse names.dmp, nodes.dmp
    nodesDB_default = join(blobtools.DATADIR, "nodesDB.txt")
    nodesDB, nodesDB_f = BtIO.parseNodesDB(nodes=nodes_f, names=names_f, nodesDB=None, nodesDBdefault=nodesDB_default)

if __name__ == '__main__':
    main()
