# job specific
# === Simplify for Standard Data processing ===
#process.ana_step.remove(process.looseTrackAna)
process.ana_step.remove(process.looseTrackAna_STD)
process.ana_step.remove(process.loosetrkEffAnalyzer)
process.ana_step.remove(process.trackAna)
process.ana_step.remove(process.refitTrackAna)
if process.trackAna.isGEN.value():
  process.schedule.remove(process.gen_step)
else:
  process.ana_step.remove(process.trkEffAnalyzer)

# === Ana Ouput Content ===
#process.trackAna_STD.histOnly = False
#process.trkEffAnalyzer.fillNtuples = True

# check spectra analyzer config
print "=== trackAna_STD cfg: ==="
print process.trackAna_STD.dumpPython()
print "=== looseTrackAna_STD cfg: ==="
print process.looseTrackAna_STD.dumpPython()
print "=== djcalo cfg: ==="
print process.djcalo.dumpPython()
print "=== djcalo_tower cfg: ==="
print process.djcalo_tower.dumpPython()

# === Addition of Customization ===
print "\nBegin My Customization\n"
# command line inputs
import sys
print "cmssw command line arguements:",sys.argv
import FWCore.ParameterSet.VarParsing as VarParsing
myopts = VarParsing.VarParsing ('standard')
# set default values
myopts.maxEvents = -1
myopts.output = 'default_output.root'
# get and parse the command line arguments
myopts.parseArguments()

# now define the customization
def mycustomise(process):
# process.myopts
  process.myopts = cms.untracked.PSet(
      wantSummary = cms.untracked.bool(True)
      )
  # process.source
  runInput=False
  for i in sys.argv:
    if ("files" in i) or("files_load" in i):
      runInput=True
      break
  if (runInput):
    process.source.fileNames = cms.untracked.vstring(myopts.files)
  # process.maxEvents
  process.maxEvents = cms.untracked.PSet(
      input = cms.untracked.int32(myopts.maxEvents)
      )
  # process.output
  try:
    if ('default_output' not in myopts.output):
      process.TFileService.fileName = cms.string(myopts.output)
  except:
    print "no output module with the given name"
  # done
  return(process)

# redefine process
process=mycustomise(process)

