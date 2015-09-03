#SVelter - Support

###Description
This folder contains two supportive files required for SVelter to run: exclude.ref.bed and CN2.ref.bed.  
`exclude.ref.bed` indicates genomic regions that would be excluded from SV analysis. Users could either use the files provided in this folder, or their own customized versions to exclude genomic regions that they are not interested in.

`CN2.ref.bed` indicates copy neutral genomic regions, where no SVs have yet been reported. Null distribution of insert size of paired end read pairs and read depth would be built based on such regions. If not provided, SVelter would randomly sample genomic regions for null model foundation.

Both files are not required but strongly recommended to be included in the *SVelter Index* step.

###Origins of each file

`exclude.ref.bed` were merged from *wgEncodeDacMapabilityConsensusExcludable.bed* and *wgEncodeDukeMapabilityRegionsExcludable.bed, both of which are publicly available* from UCSC:
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz 
```

`CN2.ref.bed` were ... ... 

*liftover* was used to transfer both files to different version accurding the reference genome used 
