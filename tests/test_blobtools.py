#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import
import unittest
from tempfile import mkdtemp
from StringIO import StringIO
from sys import path
import os
from os.path import basename, isfile, join, dirname, abspath
path.append(dirname(dirname(abspath(__file__))))
import blobtools as blobtools
import lib.BtIO as BtIO
import re
import mock

class TestBlobtoolsMain(unittest.TestCase):

    '''
    Test for cmdline behaviour
    - Look into python -munittest discover
    '''
    @mock.patch('sys.argv', ['blobtools.py', '--help'])
    #@mock.patch('sys.exit')
    @mock.patch('sys.stdout', new_callable=StringIO)
    def test_help(self, mock_stdout):
        '''
        Test : blobtools --help
        '''
        try:
            blobtools.main()
        except SystemExit:
            pass

        #print mock_stdout.getvalue()
        self.assertEqual(mock_stdout.getvalue().strip(), blobtools.__doc__.strip())


    @mock.patch('sys.stdout', new_callable=StringIO) # everything in STDOUT gets captured into mock_stdout
    def test_create(self, mock_stdout):
        with mock.patch('sys.argv', ('''
            ./blobtools create
            -i test_data/test1.nHa.fna
            -b test_data/test1.nHa.MiSeq.400.vs.assembly.mapping.bam
            -b test_data/test1.nHa.MiSeq.600.vs.assembly.mapping.bam
            -t test_data/test1.nHa.diamond.out
            -t test_data/test1.nHa.blast.out
            -o %s/foo''' % (self.tmpdir)).strip().split() ):

            try:
                blobtools.main()
            except SystemExit:
                pass
            self.assertEqual(os.listdir(self.tmpdir), ['foo.blobDB.json'])



    def setUp(self):
        self.tmpdir = mkdtemp()

    def tearDown(self):
        print("rm -rf %s" % (self.tmpdir))
        #self.tmpdir =
if __name__ == '__main__':
    unittest.main()
