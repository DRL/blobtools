import unittest
import blobtools as blobtools
import lib.BtIO as BtIO
import re

class fasta_tests(unittest.TestCase):

    '''
    Tests functions related to FASTA parsing
    >1
    AGCTAGC
    >2_spades_
    AGATGAGATGC
    >3_
    '''

    fasta_f = blobtools.TEST_DIR + "/test.fna"

    def fasta_parse(self):
        binomen_pattern = r"^([A-Z][a-z]+) ([a-z]+)$"
        self.assertEqual(re.search(binomen_pattern, 'homo sapiens'), None)

if __name__ == '__main__':
    unittest.main()
