#!/bin/bash

echo
echo "All of doctest's code is directly included in this respository, so there is usually no need to run this script,"
echo "unless you want to update doctest to the latest version. In this case, first delete the existing doctest directory, then run this script."
echo

set -ex
# In the above line, "set -e" causes this script to exit as soon as
# any line fails. "set -x" causes each line of this script to be
# printed (with a + in front) before it is executed, so if a step
# fails, you can see from the travis log what command failed.

echo Hello from install_doctest.sh

pwd

mkdir doctest
cd doctest
# The following command gets the latest version:
wget https://raw.githubusercontent.com/onqtam/doctest/master/doctest/doctest.h
# Or, you can use the following command to get a specific version:
# wget https://raw.githubusercontent.com/onqtam/doctest/2.4.3/doctest/doctest.h
cd ..

echo Done installing doctest
