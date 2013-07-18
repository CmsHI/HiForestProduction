@files = qx(ls);

@pbs = grep(/\.pbs$/,@files);

foreach $p (@pbs)
{
  chomp($p);

  system("qsub $p");
}

