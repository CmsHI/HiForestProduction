hltGetConfiguration --full --offline --data /users/frankma/devCMSSW_4_2_0/HIonV1047/V4 --type HIon --unprescale --process HLTHIon --globaltag START44_V6::All > HLT441.py

cat > tmp.py <<EOF
doL1Xml=True;
newSource=True;
if doL1Xml:
   import L1Trigger.Configuration.L1Trigger_custom
   process = L1Trigger.Configuration.L1Trigger_custom.customiseL1Menu( process )

if (newSource):
   from FWCore.ParameterSet import Mixins
   for module in process.__dict__.itervalues():
      if isinstance(module, Mixins._Parameterizable):
	 for parameter in module.__dict__.itervalues():
	    if isinstance(parameter, cms.InputTag):
	       if parameter.moduleLabel == 'source':
		  parameter.moduleLabel = 'sourcenew'
EOF

cat tmp.py >> HLT441.py
rm tmp.py
