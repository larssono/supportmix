#!/bin/bash

#This script runs SupportMix configuration tests

echo '*************** Running SupportMix Tests*************'

CONFIG_LIST="./testConfigs/incorrectFiles.cfg ./testConfigs/badLabels.cfg ./testConfigs/noAncestralFiles.cfg ./testConfigs/badRGB.cfg ./testConfigs/testConfig.cfg"

for testConfig in $CONFIG_LIST
do
   echo "************* Testing [$testConfig] ***************"
   echo "run command: ./SupportMix -t -C $testConfig"
   ./SupportMix -t -C $testConfig
   echo "***************************************************"
done
