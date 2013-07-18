#!/usr/bin/perl


$command = "gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=combinedpdf.pdf -dBATCH ";

$command .= "page_1to2.pdf ";

for ( $page = 3; $page < 66; $page++ )
{
   $command .= " page_$page.ps ";
}

system("$command");

