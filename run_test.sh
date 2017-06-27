#!/bin/bash
echo "Hello, World!" 
echo "Knowledge is power."
uname -a
sudo apt-get install build-essential
#echo $PATH
#export PATH=$PATH:/usr/local/bin/make
#echo $PATH
cd test
pwd
make verification_tests
