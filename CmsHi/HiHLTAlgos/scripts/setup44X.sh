#/bin/bash -
# hlt menu
addpkg HLTrigger/Configuration # or whatever the latest menu

# analysis examples
cvs co -d       CmsHi/HiHLTAlgos UserCode/CmsHi/HiHLTAlgos

scram build -c
scram b
