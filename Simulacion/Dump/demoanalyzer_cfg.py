import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/home/cmsusr/playground/002F62E1-B53D-E311-A49F-003048F1B950.root'
    )
)

process.demo = cms.EDAnalyzer('DemoAnalyzer'
)


process.p = cms.Path(process.demo)
