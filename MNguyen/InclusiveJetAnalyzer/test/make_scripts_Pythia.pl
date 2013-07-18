#!/usr/bin/perl



@files = `ls /pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/PatFilesWithPF_HIL1L2L3/Summer10_Pythia_Pt80/Pat_*.root`;



chomp(@files);

$i=0;
$j=0;

foreach $file (@files)
{

    #print $file."\n";


    print $outputfile;


    #if(`ls /net/hisrv0001/home/mnguyen/scratch/PatFilesWithPF_MyL2L3/AMPT_Pyquen_UnquenchedDiJet_Pt80/Pat_$i.root`)    
    #if(`ls /net/hisrv0001/home/mnguyen/scratch/PatFilesWithPF_MyL2L3/Hydjet_Pyquen_UnquenchedDiJet_Pt80/Pat_$i.root`)
    #{
    #}
    #else
    #{
	
	
    $outfile="/net/hisrv0001/home/mnguyen/scratch/CMSSW_3_9_1_patch1/src/MNguyen/InclusiveJetAnalyzer/test/exec/run_".$j++.".sh";
    open(OUT,"> $outfile") or die "$outfile\n";
    print OUT "source ~/osg_cmssw_set_basic.sh\n";
    print OUT "eval `scramv1 runtime -sh`\n";
    print OUT "cmsRun inclusiveJetAnalyzer_Pythia_cfg.py files=dcache:".$file." output=dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/InclusiveJetAnalyzer/Pythia_Pt80/jetTree_".$i.".root\n";
    
    #}
    
    $i++;
}

`chmod 777 exec/*`;
