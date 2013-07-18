import os

dcache  = "/pnfs/cmsaf.mit.edu/t2bat/cms"
scratch = "/net/hidsk0001/d00/scratch/frankma/data"
inputs_dir = '/net/hisrv0001/home/frankma/inputs_data'
samples  = {
    "MB-C10-Jun14RR-900GeV-MBskim-v1":"MinimumBias",
    "MB-C10-Jun14RR-900GeV-JMTskim-v1":"JetMETTau"
    #"MB-R10A-Jun14RR-MBskim-v0":"MinimumBias",
    #"JMT-R10A-Jun14RR-J15Uskim-v0":"JetMETTau"
    #"JMT-R10A-PR4-Jet50Uskim-v0":"JetMETTau",
    #"JMT-R10A-PR4-Jet15Uskim-v0":"JetMETTau"
    }

for sample,PD in samples.items():
  print "===== sample: ",PD,sample," =====\n"
  outdir=scratch+"/"+PD+"/"+sample
  os.system('mkdir -p '+outdir)

  filelist = []
  for line in open (inputs_dir+"/"+ sample + '.txt', 'r'):
    filelist.append(line.rstrip('\n').replace('trkAnaSkimAOD','trkhists'))

  for file in filelist:
    infile=dcache+file
    subfile = file.split('/')
    #print 'preparing to copy ' + subfile[-1] + ' from dcache'
    if not os.path.isfile(outdir + "/" + subfile[-1]):
      cmd = 'dccp ' + infile + ' ' + outdir
      print cmd
      os.system(cmd)
    else:
      print subfile[-1] + " already exists in '" + outdir + "' directory"

