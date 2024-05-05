# VertexNtuples

This is a ntuple framework for secondary vertexing studies.

### Installation with CMSSW_14_1_0_pre3
```
cmsrel CMSSW_14_1_0_pre3
cd CMSSW_14_1_0_pre3/src
cmsenv
git cms-init
git clone git@github.com:emily-tsai11/VertexNtuples.git
scram b -j
```
Or, for debugging with `gdb`, compile with
```
scram b -j USER_CXXFLAGS="-g"
```
(This is already specified by default in `VertexNtuplizer/BuildFile.xml`.)

For subsequent setup:
```
cd CMSSW_14_1_0_pre3/src/
cmsenv
```

### Running the ntuplizer locally
```
cd VertexNtuples/VertexNtuplizer/test/
cmsRun VertexNtuplizer.py inputFiles=file:/path/to/file.root
```
There is a ttbar to hadronic decays with no pileup test file here:
```
/eos/user/e/etsai/workspace/TTToHadronic_noPU_CMSSW_13_1_0/src/RecoVertex/AdaptiveVertexFinder/test/TTToHadronic_noPU_slimmed.root
```
and a ttbar to hadronic decays with 200 pileup test file here:
```
/eos/user/e/etsai/workspace/TTToHadronic_PU200_CMSSW_13_1_3/src/RecoVertex/AdaptiveVertexFinder/test/TTToHadronic_PU200_slimmed.root
```
To view the contents of these test files, use `edmDumpEventContent <filename>.root`.

### Batch submission of the ntuplizer
```
cd VertexNtuplizer/test/
crab submit crabSubmitter.py
```
