

inputDir=/home/sungho/plots/7TeV/root_files/
outputDir=./rootOutput/

if [ -d $outputDir ]
    then
    echo "$outputDir exists"
else
    mkdir $outputDir
    echo "$outputDir is created!"
fi


#spectraFile=test
#trkCorrFile1=test
#trkCorrFile2=test
#trkCorrFile3=test
#trkCorrFile4=test
#trkCorrFile5=test

#spectraFile=MB-C10-PR9-MBskim-v0_proc0628_trkAnaNoFilter_v2
spectraFile=MB-C10-PR9-JMTskim-v0_proc0628_trkAnaNoFilter_v2
trkCorrFile1=TrkHistMC_july04_qcdMBD6T0704vbinV2
trkCorrFile2=TrkHistMC_july04_qcdPt15d0704vbinV2
trkCorrFile3=TrkHistMC_july04_qcdPt15d0704vbinV2
trkCorrFile4=TrkHistMC_july04_qcdPt15d0704vbinV2
trkCorrFile5=TrkHistMC_july04_qcdPt170d0704vbinV2

#spectraFile=MB-C10-M6RR-MBskim-v0_proc0628_trkAnaNoFilter_v2_proc0628
#trkCorrFile1=TrkHistMC_june09_qcdMB
#trkCorrFile2=TrkHistMC_june09_qcdMB
#trkCorrFile3=TrkHistMC_june09_qcdMB
#trkCorrFile4=TrkHistMC_june09_qcdMB
#trkCorrFile5=TrkHistMC_june09_qcdMB


anaDir=trackAna_STD
#anaDir=trackAna
corDir=trkEffAnalyzer

isGEN=false
varBin=true

ijet=60
fjet=2000

ieta=0
feta=2.4
drawFig=true
saveFile=true

root -l CorrectTypeOneNSave.C+\(\"$inputDir\",\"$outputDir\",\"$spectraFile\",\"$trkCorrFile1\",\"$trkCorrFile2\",\"$trkCorrFile3\",\"$trkCorrFile4\",\"$trkCorrFile5\",\"$anaDir\",\"$corDir\",$isGEN\,$varBin\,$ijet\,$fjet\,$ieta\,$feta\,$drawFig\,$saveFile\)