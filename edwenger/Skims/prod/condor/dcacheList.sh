#~/bin/bash -
if [ $# -eq 0 ]; then
  echo Usage
  echo "  $0 <input_list>"
  exit 1
fi

sed -i 's;^;/pnfs/cmsaf.mit.edu/t2bat/cms;' $1
#sed -i 's;trkAnaSkimAOD;trkhists;' $1
