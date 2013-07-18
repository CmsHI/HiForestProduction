import os

dcache  = "/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/frankma"
scratch = "/net/hidsk0001/d00/scratch/frankma/data"

samples  = {
    "hcpr-spectra-v4_3-mb-all":"HCPR",
    "hcpr-spectra-v4_3-J35U-all":"HCPR",
    "hcpr-spectra-v4_3-J50U-all":"HCPR"
    }

for sample,PD in samples.items():
  print "===== sample: ",PD,sample," =====\n"
  indir=dcache+"/"+PD+"/"+sample
  outdir=scratch+"/"+PD+"/"+sample
  os.system('mkdir -p '+outdir)

  for file in os.listdir(indir):
    if "root" not in file:
      continue
    infile=indir+"/"+file
    if not os.path.isfile(outdir + "/" + file):
      cmd = 'dccp ' + infile + ' ' + outdir
      print cmd
      os.system(cmd)
    else:
      print file + " already exists in '" + outdir + "' directory"

