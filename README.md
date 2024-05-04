# VertexNTuples

Ntuple framework for vertexing studies

### Installation with CMSSW_14_1_0_pre3
```
cmsrel CMSSW_14_1_0_pre3
cd CMSSW_14_1_0_pre3/src
cmsenv
<install this ntuplizer>
```

### Compile with debug option
```
scram b -j USER_CXXFLAGS="-g"
```

### Running the ntuplizer
```
cd VertexNTuplizer/test/
<some python config>
```
