#Support files for human genome hg19 and hg38

Under this folder, we provide supportive files for running SVelter with human reference hg19 and hg38.

Each folder (except for *ref-index*) contains following files, and the path would be required by svelter Setup to properly set the working directory.

_Exclude.*.bed_

_CN2.*.bed_

_Segdup.*.bed_

_SVelter.*.r_


folder *ref-index* contains pre-indexed files for human reference hg19 and hg38. Feeding the path would accelerate svelter Setup.

####example of setting up working directory with all supportive files feeded:
```
svelter.py Setup --workdir ./workdir/ --reference ref.fa --support ../Support/hg19/ --ref-index ../Support/ref-index/hg19/
```
