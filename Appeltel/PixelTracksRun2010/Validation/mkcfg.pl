#!/usr/bin/perl

open( TMPL,"validation_CFG_template") or die $!;
@cfg = <TMPL>;
close TMPL;

open( PTMPL,"validation_PBS_template") or die $!;
@pbs = <PTMPL>;
close PTMPL;

@centstrings = ("[2,3]","[4,5]","[6,7]","[8,9,10,11]","[12,13,14,15]",
               "[16,17,18,19]","[20,21,22,23]","[24,25,26,27]",
               "[28,29,30,31]" );

@etas = ( 0.5, 0.8, 1.0, 1.5, 2.0, 2.4 );

foreach $eta (@etas)
{
  $cbin = 0;
  foreach $cent (@centstrings)
  {
     $cbin++;
     $cfgname = "val" . "_eta$eta" . "_cbin$cbin" . ".py";
     open( CFG, ">$cfgname");
     foreach $line (@cfg)
     {
       $t = $line;
       $t =~ s/#_CENT_#/$cent/g;
       $t =~ s/#_ETA_#/$eta/g;
       $t =~ s/#_CBIN_#/$cbin/g;
       print CFG $t;
     }
     close CFG;

     open ( PBS, ">val" . "_eta$eta" . "_cbin$cbin" . ".pbs" );
     foreach $line (@pbs)
     {
       $t = $line;
       $t =~ s/#_CFGNAME_#/$cfgname/g;
       $t =~ s/#_ETA_#/$eta/g;
       $t =~ s/#_CBIN_#/$cbin/g;
       print PBS $t;
     }
     close PBS;

  }
}


