inputDir=/Users/andreyoon/Research/pbpbSpectra/root_files/mc/
outputDir=./rootOutput_trk/may122011_invstg/

if [ -d $outputDir ]
    then
    echo "$outputDir exists"
else
    mkdir $outputDir
    echo "$outputDir is created!"
fi

trkCorrFile1=trkhist_may112011_v1_hydjetBass_dijet170_finec     

anaDir=hitrackAna
corDir=hitrkEffAnalyzer
#corDir=higloosetrkEffAnalyzer
#corDir=hihightrkEffAnalyzer


function runEffMacro()
{
    root -l -q GetEffAndFakeNSave.C++\(\"$inputDir\",\"$outputDir\",\"$trkCorrFile1\",\"$corDir\",$ijet\,$fjet\,$ieta\,$feta\,$ipt,$drawFig\,$saveFile\,$var\,$rebF\,$1,$2,$effonly\)
}


effonly=false

ijet=0
fjet=2000


ieta=0.0
feta=1.0

ipt=0.0
#ipt=2.0

rebF=2

drawFig=true
saveFile=true

var=1

cbini=0
cbinf=1

runEffMacro 0 1
runEffMacro 2 3
runEffMacro 4 11
runEffMacro 12 19
#runEffMacro 20 35

