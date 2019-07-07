# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: analyze --data --conditions 101X_dataRun2_Prompt_v11 -s RAW2DIGI,L1Reco,RECO --geometry DB:Extended --era Run2_2018 --eventcontent FEVTDEBUG --customise Debug/HcalDebug/customize.compare_raw_reemul_tp --customise Debug/HcalDebug/customize.use_data_reemul_tp --customise_commands cms.untracked.vstring('keep *,keep FEDRawData_*_*_*') --filein=filelist:filelist.txt
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('reco1',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(150)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/data/users/snabili/CMSSW_10_3_2/src/analyze_RAW2DIGI_L1Reco_RECO_n400.root'),
    fileNames = cms.untracked.vstring('file:/data/users/snabili/CMSSW_10_3_2/src/filter.root'),
    #fileNames = cms.untracked.vstring('file:/data/users/snabili/CMSSW_10_3_2/src/FILTER_zerobias.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('analyze nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        #filterName = cms.untracked.string('filter_step')
    ),
    fileName = cms.untracked.string('analyze_reco2.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
    #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('filter_step') )
)

# Additional output definition





# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v11', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
#process.L1Reco_step = cms.Path(process.L1Reco)
#process.reconstruction_step = cms.Path(process.reconstruction)
#process.filter_step = cms.Path(process.HBHENoiseFilter)
#process.p_step = cms.Path(
  #process.calolocalreco *
  #process.caloTowersRec *
  #process.hcalnoise *
  #process.particleFlowCluster 
#)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.FEVTDEBUGoutput_step)
#process.schedule = cms.Schedule(process.raw2digi_step,process.p_step,process.endjob_step,process.FEVTDEBUGoutput_step)
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from Debug.HcalDebug.customize
from Debug.HcalDebug.customize import compare_raw_reemul_tp,use_data_reemul_tp,compare_raw_reco_sev9 

#call to customisation function compare_raw_reemul_tp imported from Debug.HcalDebug.customize
process = compare_raw_reemul_tp(process)

#call to customisation function use_data_reemul_tp imported from Debug.HcalDebug.customize
process = use_data_reemul_tp(process)

#call to customisation function compare_reemul_reco_sev9 imported from Debug.HcalDebug.customize
process = compare_raw_reco_sev9(process)

# End of customisation functions

# Customisation from command line
process.TFileService.fileName=cms.string('reco_filt_test1.root')
cms.untracked.vstring('keep recoPFRecHits_particleFlow_*_*_*')
cms.untracked.vstring('keep *,keep FEDRawData_*_*_*')
#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
