import FWCore.ParameterSet.Config as cms

import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("acausal")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#Data files
#files = FileUtils.loadListFromFile("data/CMS_Run2011A_DoubleElectron_AOD_12Oct2013-v1_20000_file_index.txt")
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2011A_DoubleElectron_AOD_12Oct2013-v1_20001_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2011B_DoubleElectron_AOD_12Oct2013-v1_00000_file_index.txt"))
#files.extend(FileUtils.loadListFromFile("data/CMS_Run2011B_DoubleElectron_AOD_12Oct2013-v1_20000_file_index.txt")) 

#Background files example
#files = FileUtils.loadListFromFile("data/CMS_MonteCarlo2011_Summer11LegDR_DYJetsToLL_M-50_7TeV-madgraph-pythia6-tauola_AODSIM_PU_S13_START53_LV6-v1_00000_file_index.txt")

#process.source = cms.Source("PoolSource",
#    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring(
#       *files
#    )
#)

#Signal
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:data/DoubleEle.root'
    )
)

process.acausal = cms.EDAnalyzer('AcausalAnalyzer', isData = cms.bool(True)
)


#needed to get the actual prescale values used from the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')
process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'


#configure the analyzer
#inspired by https://github.com/cms-sw/cmssw/blob/CMSSW_5_3_X/HLTrigger/HLTfilters/interface/HLTHighLevel.h
process.gettriggerinfo = cms.EDAnalyzer('AcausalAnalyzer')
#                              processName = cms.string("HLT"),
#                              triggerPatterns = cms.vstring("HLT_DoubleEle45_CaloIdL_v*"), #if left empty, all triggers will run        
#                              triggerResults = cms.InputTag("TriggerResults","","HLT"),
#                              triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")                             
#                              )

process.p = cms.Path(process.acausal)
