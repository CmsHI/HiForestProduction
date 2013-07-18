import os

scratch="/net/hidsk0001/d00/scratch/frankma/data"
jobTag="proc0210_AllCent"

# ================== 39X ===============
hisamples  = {
    "hcpr-spectra-v4_3-mb-all":"HCPR",
    "hcpr-spectra-v4_3-J35U-all":"HCPR",
    "hcpr-spectra-v4_3-J50U-all":"HCPR"
    }

# ================== 362 ===============
ppsamples  = {
    # 900GeV
    #"MB-C10-Jun14RR-900GeV-MBskim-v1":"MinimumBias",
    #"MB-C10-Jun14RR-900GeV-JMTskim-v1":"JetMETTau",
    # 7TeV
    "MB-C10-M6RR-JMTskim-v0":"MinimumBias",
    "MB-C10-PR9-JMTskim-v0":"MinimumBias",
    "JMT-R10A-PR1-JMTskim-v0":"JetMETTau",
    "JM-R10A-PR4-Jet50Uskim-v0":"JetMETTau",
    "JM-R10A-PR4-Jet15Uskim-v0":"JetMETTau"
    }

inFolder="."
outFolder="mergeAll"
for sample,PD in hisamples.items():
  #indir=scratch+"/"+PD+"/"+sample+"_"+jobTag
  indir=scratch+"/"+PD+"/"+sample
  cmd="./mergeJob.sh " + indir+"/"+inFolder + " trkhists " + indir+"/"+outFolder
  print cmd
  os.system(cmd)
