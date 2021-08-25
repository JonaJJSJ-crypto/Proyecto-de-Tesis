import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter("Pythia8HadronizerFilter",
   maxEventsToPrint = cms.untracked.int32(0),
   pythiaPylistVerbosity = cms.untracked.int32(1),
   filterEfficiency = cms.untracked.double(1.0),
   pythiaHepMCVerbosity = cms.untracked.bool(True),
   SLHAFileForPythia8 = cms.string('SimLW/Sim/param_card.slha'),
   comEnergy = cms.double(8000.0),
   UseExternalGenerators = cms.untracked.bool(True),
   PythiaParameters = cms.PSet(
       processParameters = cms.vstring(
           'Tune:pp 5',
           'PDF:pSet = 5',
           '111:mayDecay = off'
       ),
       parameterSets = cms.vstring('processParameters')
   )
)
