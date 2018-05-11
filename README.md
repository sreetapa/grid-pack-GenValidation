# grid-pack-GenValidation
cerate cmsenv with specific release and put the directory in CMSSW/src

do "scram b"

change GEN-SIM input file name in test/runAnalyzer.py  

do "cmsRun runAnalyzer.py"

it will create a root file with basic GenJet kinematic distribution
