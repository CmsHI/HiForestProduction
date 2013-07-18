menu=$1
cfg=$2
cmd="hltGetConfiguration --full --offline --data $menu --type HIon --input root://eoscms//eos/cms/store/hidata/HIRun2011/HIMinBiasUPC/RAW/v1/000/182/296/F6A6C5C7-4615-E111-8E79-003048F11DE2.root --process HLT2 > $cfg"
echo $cmd
eval $cmd
sed -i 's/"source"/"rawDataRepacker"/g' $cfg
sed -i "s/'source'/'rawDataRepacker'/g" $cfg
