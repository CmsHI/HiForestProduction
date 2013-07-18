#!/usr/bin/perl



@files = `ls /pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/HICorePhysics/Jet50-PromptReco-Runs_151350*/*/RECOPAT*.root`;



chomp(@files);



$i=0;
$j=0;

foreach $file (@files)
{

    #print $file."\n";
    @runfile = split(/\//,$file);

    print $runfile[9]."\n";
    print $runfile[11]."\n";

    @runextra = split(/Jet50-PromptReco-/,$runfile[9]);
    @run = split(/_RECOPAT/,$runextra[1]);

    $outputfile = "jetTree_MoreAlgos_".$run[0]."_".$runfile[11];


    if(`ls /pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/InclusiveJetAnalyzer/HIData_Jet50U/$outputfile`)
	{
	}
    else
    {
	
	
    $outfile="/net/hisrv0001/home/mnguyen/scratch/CMSSW_3_9_1_patch1/src/MNguyen/InclusiveJetAnalyzer/test/exec/run_".$j++.".sh";
    open(OUT,"> $outfile") or die "$outfile\n";
    print OUT "source ~/osg_cmssw_set_basic.sh\n";
    print OUT "eval `scramv1 runtime -sh`\n";
    print OUT "cmsRun inclusiveJetAnalyzer_Data_cfg.py files=dcache:".$file." output=dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/mnguyen/InclusiveJetAnalyzer/HIData_Jet50U/".$outputfile."\n";
    
    }
    
    $i++;
}

`chmod 777 exec/*`;
