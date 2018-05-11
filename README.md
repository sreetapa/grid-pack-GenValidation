# grid-pack-GenValidation
cerate cmsenv with specific CMSSW release and go to src area

do "git clone git@github.com:panwarlsweet/grid-pack-GenValidation.git GenValidation/NtupleGenJet"

do "scram b" in src area 

change GEN-SIM input file name in test/runAnalyzer.py  

do "cmsRun runAnalyzer.py"

it will create a root file with basic GenJet kinematic distribution
