import os

scratch="/net/hidsk0001/d00/scratch/frankma"
jobTag="proc0111_JetVtx"

# ================== 362 ===============
samples  = {
    #"MB-C10-M6RR-JMTskim-v0":"MinimumBias",
    #"MB-C10-PR9-JMTskim-v0":"MinimumBias",
    #"JMT-R10A-PR1-JMTskim-v0":"JetMETTau",
    #"JM-R10A-PR4-Jet50Uskim-v0":"JetMETTau",
    "JM-R10A-PR4-Jet15Uskim-v0":"JetMETTau"
    }

for sample,PD in samples.items():
  outputDir=scratch+"/data/"+PD+"/"+sample+"_"+jobTag
  print outputDir
  os.chdir(outputDir)
  os.system('echo "  # of files:" `ls trkhists*.root | wc -l`')
  cmd="../../scripts/cleanJob.sh"
  print " run: ",cmd
  os.system(cmd)

