import FWCore.ParameterSet.Config as cms

process = cms.Process("Validation")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
            'file:../../../../../2016_LHEs/HIG-RunIISummer15wmLHEGS-01479.root'
                )
                            )


ntuple_genjet = cms.PSet(
    NtupleName = cms.string('NtupleGenJet'),
    GenJets = cms.InputTag('ak4GenJets')
)

process.demo = cms.EDAnalyzer(
    "NtupleGenJet",
    Ntuples = cms.VPSet(
        ntuple_genjet,
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("out_tree_susy_M200_2016.root"
))
process.p = cms.Path(process.demo)
process.MessageLogger.cerr.FwkReport.reportEvery = 100


