#!/usr/bin/perl



#@files = `ls /pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/PatFilesWithPF_HIL1L2L3/Hydjet_Pyquen_DiJet_Pt80MC_38Y_V12-v2/Pat_*.root`;
#@files = `ls /pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/PatFilesWithPF_HIL1L2L3/AMPT_Pyquen_DiJet_Pt80_MC_38Y_V12-v1/Pat_*.root`;
@files = `ls /pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/PatFilesWithPF_HIL1L2L3/Hydjet_Pyquen_UnquenchedDiJet_Pt80/Pat_FJ2Jets*.root`;


chomp(@files);

$i=0;
$j=0;

foreach $file (@files)
{

    #print $file."\n";


    print $outputfile;


    #if(`ls /pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/InclusiveJetAnalyzer/Hydjet_Pyquen_Pt80MC_38Y_V12-v2/jetTree_$i.root`)
    #if(`ls /pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/InclusiveJetAnalyzer/AMPT_Pyquen_Pt80_MC_38Y_V12-v1/jetTree_$i.root`)
    if(`ls /pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/InclusiveJetAnalyzer/Hydjet_Pyquen_Unquenched_Pt80/jetTree_MoreAlgos_$i.root`)
    {
    }
    else
    {
	
	
    $outfile="/net/hisrv0001/home/mnguyen/scratch/CMSSW_3_9_1_patch1/src/MNguyen/InclusiveJetAnalyzer/test/exec/run_".$j++.".sh";
    open(OUT,"> $outfile") or die "$outfile\n";
    print OUT "source ~/osg_cmssw_set_basic.sh\n";
    print OUT "eval `scramv1 runtime -sh`\n";
    print OUT "cmsRun inclusiveJetAnalyzer_Embedded_cfg.py files=dcache:".$file." output=dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/InclusiveJetAnalyzer/Hydjet_Pyquen_Unquenched_Pt80/jetTree_MoreAlgos_".$i.".root\n";
    
    }
    
    $i++;
}

`chmod 777 exec/*`;
