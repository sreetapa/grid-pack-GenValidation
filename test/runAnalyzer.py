import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the sourcefile you want to use
                                fileNames = cms.untracked.vstring(
#'/store/user/gkole/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190817_164543/0001/B2G-RunIIFall18GS-00005_3_1000.root',
#'/store/user/gkole/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190817_164543/0000/B2G-RunIIFall18GS-00005_3_999.root ',
#'/store/user/gkole/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190817_164543/0000/B2G-RunIIFall18GS-00005_3_998.root',
#'/store/user/gkole/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190817_164543/0000/B2G-RunIIFall18GS-00005_3_997.root',
#'/store/user/gkole/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190817_164543/0000/B2G-RunIIFall18GS-00005_3_996.root',
#'/store/user/gkole/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190817_164543/0000/B2G-RunIIFall18GS-00005_3_995.root'
#'/store/user/spoddar/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190827_123049/0000/B2G-RunIIFall18GS-00005_3_999.root'
#'/store/user/spoddar/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190828_121700/0000/B2G-RunIIFall18GS-00005_3_175.root',
#'/store/user/spoddar/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190828_121700/0000/B2G-RunIIFall18GS-00005_3_534.root',
#'/store/user/spoddar/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190828_121700/0000/B2G-RunIIFall18GS-00005_3_533.root',
#'/store/user/spoddar/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190828_121700/0001/B2G-RunIIFall18GS-00005_3_1000.root'
#"file:/afs/cern.ch/work/s/spoddar/private/CMSSW_10_2_3/src/crab_MC_Generation/B2G-RunIIFall18GS-00005_3_mass_700_cfg.root"
#'/store/user/spoddar/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190903_061149/0000/B2G-RunIIFall18GS-00005_3_8.root'
'file:/afs/cern.ch/work/s/spoddar/private/CMSSW_10_2_3/src/crab_MC_Generation/B2G-RunIIFall18GS-00005_3.root'

#'file:/afs/cern.ch/work/s/spoddar/private/CMSSW_10_2_6/src/GenValidation/NtupleGenJet/test/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190827_121700/0000/B2G-RunIIFall18GS-00005_3_997.root'

#"file:/afs/cern.ch/work/s/spoddar/private/CMSSW_10_2_6/src/GenValidation/NtupleGenJet/test/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190828_121700/0000/B2G-RunIIFall18GS-00005_3_996.root'
#'file:/afs/cern.ch/work/s/spoddar/private/CMSSW_10_2_6/src/GenValidation/NtupleGenJet/test/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190828_121700/0000/B2G-RunIIFall18GS-00005_3_995.root'
#'file:/afs/cern.ch/work/s/spoddar/private/CMSSW_10_2_6/src/GenValidation/NtupleGenJet/test/TprimeTPrimeTotGammatGluon/TprimeTPrimeTotGammatGluon_cfg_GEN-SIM/190828_121700/0000/B2G-RunIIFall18GS-00005_3_994.root'






        #'file:/afs/cern.ch/work/l/lata/grid_pack/genproductions/bin/JHUGen/production_validation_HIG/susy_gridpack_validation/2016_LHEs/HIG-RunIISummer15wmLHEGS-01479.root'
        
#'file:/afs/cern.ch/user/g/gkole/work/public/HIG-RunIISummer15wmLHEGS-01479.root' # for bbHToTauTau 
#'file:/eos/user/g/gkole/PostDoc/temp/B2G-RunIIFall18GS-00005_3_nEvents500.root' # for tprimeToTGamma

                )
                            )



ntuple_genHiggs = cms.PSet(
     NtupleName = cms.string('NtupleGenJet'),
     genParticles = cms.InputTag('genParticles'),
)
process.demo = cms.EDAnalyzer(
    "NtupleGenJet",
    Ntuples = cms.VPSet(
	ntuple_genHiggs,
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out_tree_TPrime_M700_2018.root"
))
process.p = cms.Path(process.demo)
process.MessageLogger.cerr.FwkReport.reportEvery = 100


