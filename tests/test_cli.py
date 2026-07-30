#!/usr/bin/env python3

"""Tests for the xbooz_xform command-line driver (booz_xform._cli)."""

import os
import shutil
import tempfile
import unittest

from booz_xform import _cli

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_files')


class ParseInputFileTest(unittest.TestCase):
    def _parse(self, contents):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(tmpdir, 'in_booz.test')
            with open(filename, 'w') as f:
                f.write(contents)
            return _cli._parse_input_file(filename)

    def test_with_compute_surfs(self):
        mboz, nboz, extension, compute_surfs = self._parse(
            "32 16\nli383_1.4m\n5 10 15\n")
        self.assertEqual(mboz, 32)
        self.assertEqual(nboz, 16)
        self.assertEqual(extension, 'li383_1.4m')
        self.assertEqual(compute_surfs, [5, 10, 15])

    def test_without_compute_surfs(self):
        mboz, nboz, extension, compute_surfs = self._parse("32 16\nli383_1.4m\n")
        self.assertEqual((mboz, nboz, extension), (32, 16, 'li383_1.4m'))
        self.assertEqual(compute_surfs, [])

    def test_too_few_entries(self):
        with self.assertRaises(RuntimeError):
            self._parse("32 16\n")


class RunDriverTest(unittest.TestCase):
    def test_transformation_from_the_command_line(self):
        """Run the driver end to end, as the xbooz_xform script does.

        The driver reads wout_<extension>.nc and writes boozmn_<extension>.nc,
        both relative to the working directory.
        """
        configuration = 'circular_tokamak'
        original_dir = os.getcwd()
        with tempfile.TemporaryDirectory() as tmpdir:
            shutil.copy(
                os.path.join(TEST_DIR, 'wout_' + configuration + '.nc'), tmpdir)
            with open(os.path.join(tmpdir, 'in_booz'), 'w') as f:
                f.write("16 8\n" + configuration + "\n0 1\n")
            try:
                os.chdir(tmpdir)
                self.assertEqual(_cli.main(['in_booz']), 0)
                self.assertTrue(
                    os.path.isfile('boozmn_' + configuration + '.nc'))
            finally:
                os.chdir(original_dir)

    def test_usage_message(self):
        """No arguments prints usage and exits successfully, as in driver.cpp."""
        self.assertEqual(_cli.main([]), 0)


if __name__ == "__main__":
    unittest.main()
