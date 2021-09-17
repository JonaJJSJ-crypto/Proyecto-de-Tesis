# Auto generated configuration file
# using: 
# Revision: 1.381.2.28 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: SimLW/Sim/python/myHadronizer.py --step GEN,SIM --datatier GEN-SIM --conditions START53_V27::All --beamspot Realistic8TeVCollision --customise Configuration/DataProcessing/Utils.addMonitoring --fileout gensimLW.root --eventcontent RAWSIM --no_exec --filein file:LWSM200DnR.lhe --filetype LHE --number=10 --python_filename gensimLW.py --mc
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring('file:LWSM200DnR.lhe'),
    inputCommands = cms.untracked.vstring('keep *', 
        'drop LHEXMLStringProduct_*_*_*'),
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.28 $'),
    annotation = cms.untracked.string('SimLW/Sim/python/myHadronizer.py nevts:10'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('gensimLW.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V27::All', '')
#process.GlobalTag.connect = cms.string('sqlite_file:/opt/cms-opendata-conddb/START53_V27_MC_stripped.db')
#process.GlobalTag.globaltag = "START53_V27::All"


process.generator = cms.EDFilter("Pythia8HadronizerFilter",
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    comEnergy = cms.double(8000.0),
    UseExternalGenerators = cms.untracked.bool(True),
    PythiaParameters = cms.PSet(
        processParameters = cms.vstring('Tune:pp 5', 
            'PDF:pSet = 5'),
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
        
    ),
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
