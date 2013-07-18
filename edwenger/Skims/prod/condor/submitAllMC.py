import os
os.system("./create_custom_python_cfg.sh")

scratch="/net/hidsk0001/d00/scratch/frankma"
jobTag="proc1022_final_LIGHT_v3"

# ================== 362 ===============
samples  = {
    "su10-qcdflat-startup36v9":"QCD"
    #"su10-qcd80-startup36v9_f500":"QCD"
    #"su10-qcd80to120":"QCD"
    }

for sample,PD in samples.items():
  inputList="~/inputs_mc/"+sample+".txt"
  outputDir=scratch+"/mc/"+PD+"/"+sample+"_"+jobTag
  cmd="./submitJob.sh "+inputList+" "+outputDir
  print cmd
  os.system(cmd)
