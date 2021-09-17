import FWCore.ParameterSet.Config as cms

generator = cms.EDFilter("Pythia8HadronizerFilter",
   maxEventsToPrint = cms.untracked.int32(0),
   pythiaPylistVerbosity = cms.untracked.int32(1),
   filterEfficiency = cms.untracked.double(1.0),
   pythiaHepMCVerbosity = cms.untracked.bool(True),
   comEnergy = cms.double(8000.0),
   UseExternalGenerators = cms.untracked.bool(True),
   PythiaParameters = cms.PSet(
       processParameters = cms.vstring(
           'Tune:pp 5',
           'PDF:pSet = 5',
       ),
       pythiaMyParameters = cms.vstring(
         '556:new = lwe- lwe+ 2 -3 0 200.0 0.0 200.0 200.0 2.70765e-02',
         '556:isResonance=off',
         '556:isVisible=off',
         '556:addChannel= 1 1.0 100 23 11',
         '23:isResonance=off',
	 '23:oneChannel= 1   0.1540492    0        1	   -1',
         '23:addChannel= 1   0.1194935    0        2	   -2',
         '23:addChannel= 1   0.1540386    0        3	   -3',
         '23:addChannel= 1   0.1193325    0        4	   -4',
         '23:addChannel= 1   0.1523269    0        5	   -5'
        ),

	    parameterSets = cms.vstring('processParameters',
                                    'pythiaMyParameters')
   )
)
