#!/usr/bin/env python3

import unittest
import os
from booz_xform import Booz_xform

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_files')

class MainTest(unittest.TestCase):
    def test_jlist_edit(self):
        print("dir: ", dir())
        b = Booz_xform()
        b.read_wout(os.path.join(TEST_DIR, 'wout_li383_1.4m.nc'))
        b.jlist = [10, 15]
        self.assertEqual(b.jlist, [10, 15])

if __name__ == '__main__':
    unittest.main()
