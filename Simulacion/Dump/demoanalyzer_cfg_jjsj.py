import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'root://eospublic.cern.ch//eos/opendata/cms/Run2012C/DoubleElectron/AOD/22Jan2013-v1/20000/0032A48E-EA67-E211-AC23-0026189438BD.root'
    )
)

process.demo = cms.EDAnalyzer('DemoAnalyzer'
)


process.p = cms.Path(process.demo)
