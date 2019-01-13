# grid-pack-GenValidation
cerate cmsenv with specific CMSSW release and go to src area

Instructions for installation:
 ```
 cmsrel CMSSW_10_2_6
 cd CMSSW_10_2_6/src
 cmsenv
 git clone -b master_mine git@github.com:gourangakole/grid-pack-GenValidation.git GenValidation/NtupleGenJet
 scram b
 ```

change GEN-SIM input file name in test/runAnalyzer.py  

do "cmsRun runAnalyzer.py"

it will create a root file with basic GenJet kinematic distribution
