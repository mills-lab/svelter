#Compare SV detecting Results from different algorithms:

Four different algortihms: *SVelter*, *Delly*, *Lumpy* and *Pindel* were applied on simulated data and compared;
we used several **shell**, **python** and **R** commands to accomplish the whole process. 

##scripts used:
`SV.Output.Process.py `

`Produce.Pseudo.ROC.stats.py`

###Usage 
SV.Output.Process.py  [options] [parameter]
####Options:
```
vcf-to-bed:  extract simple SVs from vcf files and output in separate bed files
bedpe-to-bed:    extract simple SVs from bedpe files and output in separate bed files
Mappable-Control:    remove SVs located outside mappable regions
Size-Control:    filter out SVs of size outside defined range
TRA-Control: remove SVs overlap with defined SVs
```
####Required Parameters:
######Parameters for vcf-to-bed:
```
--input: input file
``` 
######Parameters for bedpe-to-bed:
```
--input: input file
--reference: reference.genome.fa
```
######Parameters for Mappable-Control:
```
--input: input file
--ref-prefix: reference.genome.fa
```
######Parameters for Size-Control:
```
--input: input.bed
--min-size: reference.Mappable. describing mappable regions
--max-size: reference.Mappable. describing mappable regions
```
######Parameters for TRA-Control:
```
--input: input.bed
--TRA-rec: TRA information kept in .rec files
```

##Example Workflow:
>####Step1. Remove low quality calls from vcf / bedpe files
```
grep -v LowQual Input.vcf>Input_QC.vcf
```
####Step2. Extract simple SVs from vcf/bedpe file to write in bed format.
```
SV.Output.Process.py vcf --input Input_QC.vcf
SV.Output.Process.py bedpe --input Input_QC.vcf --reference genome.fa
```
####Step3. Extract only SVs within mappable gemoic regions:
```
SV.Output.Process.py Mappable-Control --input Input_QC.SV.bed --ref-prefix genome.Mappable
```
####Step4. Remove any predictions that interset with preset Translocations
as other algorithms we compapre SVelter to tend to misinteprete translocation as deletion+duplication. 
```
SV.Output.Process.py TRA-Control --TRA-rec Simulate.TRA.rec --input Input_QC.SV.Mappable.bed
```
####Step5. Keep SVs within expected size range
```
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD10_sorted.DUP.Mappable.TRAFree.bed
```
####Step6. Compare predicted SVs against simulated set to decide sensitivity and specificity of each algorithm:
```
Produce.Pseudo.ROC.stats.py --path_ref folder/contains/simulatd.sv.bed --path_in folder/contains/predicted.sv.bed --appdix .Mappable.TRAFree.min100.max1000000000.bed
```
###Step7. Visualize
```
```


###Detailed Code Used:
>!
#####To process **Simulated Simple Homo** events:
######SVelter:
```
change.vcf.name.py -p ./

grep -v LowQual homo_RD10_sorted.vcf>SVelter_QC_homo_RD10_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_homo_RD10_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD10_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD10_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD10_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD10_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD10_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD10_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD10_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD10_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD10_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD10_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD10_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD10_sorted.DUP.Mappable.TRAFree.bed

grep -v LowQual homo_RD20_sorted.vcf > SVelter_QC_homo_RD20_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_homo_RD20_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD20_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD20_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD20_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD20_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD20_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD20_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD20_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD20_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD20_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD20_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD20_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD20_sorted.DUP.Mappable.TRAFree.bed

grep -v LowQual homo_RD30_sorted.vcf > SVelter_QC_homo_RD30_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_homo_RD30_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD30_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD30_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD30_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD30_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD30_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD30_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD30_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD30_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD30_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD30_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD30_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD30_sorted.DUP.Mappable.TRAFree.bed

grep -v LowQual homo_RD40_sorted.vcf > SVelter_QC_homo_RD40_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_homo_RD40_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD40_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD40_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD40_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD40_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD40_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD40_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD40_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD40_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD40_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD40_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD40_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD40_sorted.DUP.Mappable.TRAFree.bed

grep -v LowQual homo_RD50_sorted.vcf > SVelter_QC_homo_RD50_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_homo_RD50_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD50_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD50_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD50_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_homo_RD50_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD50_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD50_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD50_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input SVelter_QC_homo_RD50_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD50_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD50_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD50_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_homo_RD50_sorted.DUP.Mappable.TRAFree.bed
```
######Delly
```
change.vcf.name.py -p ./

grep -v LowQual homo_RD10_DEL.vcf > Delly_QC_homo_RD10_DEL.vcf
grep -v LowQual homo_RD10_DUP.vcf > Delly_QC_homo_RD10_DUP.vcf
grep -v LowQual homo_RD10_INV.vcf > Delly_QC_homo_RD10_INV.vcf
grep -v LowQual homo_RD10_TRA.vcf > Delly_QC_homo_RD10_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_homo_RD10_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD10_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD10_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD10_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD10_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD10_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD10_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD10_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD10_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD10_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD10_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD10_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD10_INV.INV.Mappable.TRAFree.bed

grep -v LowQual homo_RD20_DEL.vcf > Delly_QC_homo_RD20_DEL.vcf
grep -v LowQual homo_RD20_DUP.vcf > Delly_QC_homo_RD20_DUP.vcf
grep -v LowQual homo_RD20_INV.vcf > Delly_QC_homo_RD20_INV.vcf
grep -v LowQual homo_RD20_TRA.vcf > Delly_QC_homo_RD20_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_homo_RD20_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD20_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD20_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD20_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD20_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD20_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD20_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD20_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD20_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD20_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD20_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD20_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD20_INV.INV.Mappable.TRAFree.bed

grep -v LowQual homo_RD30_DEL.vcf > Delly_QC_homo_RD30_DEL.vcf
grep -v LowQual homo_RD30_DUP.vcf > Delly_QC_homo_RD30_DUP.vcf
grep -v LowQual homo_RD30_INV.vcf > Delly_QC_homo_RD30_INV.vcf
grep -v LowQual homo_RD30_TRA.vcf > Delly_QC_homo_RD30_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_homo_RD30_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD30_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD30_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD30_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD30_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD30_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD30_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD30_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD30_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD30_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD30_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD30_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD30_INV.INV.Mappable.TRAFree.bed

grep -v LowQual homo_RD40_DEL.vcf > Delly_QC_homo_RD40_DEL.vcf
grep -v LowQual homo_RD40_DUP.vcf > Delly_QC_homo_RD40_DUP.vcf
grep -v LowQual homo_RD40_INV.vcf > Delly_QC_homo_RD40_INV.vcf
grep -v LowQual homo_RD40_TRA.vcf > Delly_QC_homo_RD40_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_homo_RD40_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD40_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD40_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD40_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD40_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD40_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD40_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD40_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD40_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD40_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD40_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD40_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD40_INV.INV.Mappable.TRAFree.bed

grep -v LowQual homo_RD50_DEL.vcf > Delly_QC_homo_RD50_DEL.vcf
grep -v LowQual homo_RD50_DUP.vcf > Delly_QC_homo_RD50_DUP.vcf
grep -v LowQual homo_RD50_INV.vcf > Delly_QC_homo_RD50_INV.vcf
grep -v LowQual homo_RD50_TRA.vcf > Delly_QC_homo_RD50_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_homo_RD50_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD50_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD50_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_homo_RD50_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD50_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD50_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_homo_RD50_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD50_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD50_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Delly_QC_homo_RD50_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD50_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD50_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_homo_RD50_INV.INV.Mappable.TRAFree.bed

ln -s Delly_QC_homo_RD10_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_homo_RD10_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Delly_QC_homo_RD20_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_homo_RD20_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Delly_QC_homo_RD30_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_homo_RD30_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Delly_QC_homo_RD40_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_homo_RD40_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Delly_QC_homo_RD50_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_homo_RD50_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
```
######Lumpy
```
change.bedpe.name.py -p ./

mv homo_RD10_pesr.bedpe Lumpy_QC_RD10_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD10_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD10_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD10_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD10_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD10_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD10_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD10_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD10_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD10_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD10_pesr.INV.Mappable.TRAFree.bed

mv homo_RD20_pesr.bedpe Lumpy_QC_RD20_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD20_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD20_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD20_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD20_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD20_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD20_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD20_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD20_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD20_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD20_pesr.INV.Mappable.TRAFree.bed

mv homo_RD30_pesr.bedpe Lumpy_QC_RD30_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD30_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta 
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD30_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD30_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD30_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD30_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD30_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD30_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD30_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD30_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD30_pesr.INV.Mappable.TRAFree.bed

mv homo_RD40_pesr.bedpe Lumpy_QC_RD40_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD40_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta 
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD40_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD40_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD40_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD40_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD40_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD40_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD40_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD40_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD40_pesr.INV.Mappable.TRAFree.bed

mv homo_RD50_pesr.bedpe Lumpy_QC_RD50_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD50_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta 
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD50_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD50_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD50_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD50_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD50_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Lumpy_QC_RD50_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD50_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD50_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD50_pesr.INV.Mappable.TRAFree.bed

ln -s Lumpy_QC_RD10_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD10_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Lumpy_QC_RD20_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD20_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Lumpy_QC_RD30_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD30_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Lumpy_QC_RD40_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD40_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Lumpy_QC_RD50_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD50_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
```
######Pindel
```
change.vcf.name.py -p ./

vcf.size.filter.py -i Pindel_homo_RD_10.vcf --size 100
grep -v LowQual Pindel_homo_RD_10.LargerThan100.vcf > Pindel_QC_homo_RD10_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_homo_RD10_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD10_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD10_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD10_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD10_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD10_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD10_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD10_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD10_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD10_LargerThan100.INV.Mappable.TRAFree.bed

vcf.size.filter.py -i Pindel_homo_RD_20.vcf --size 100
grep -v LowQual Pindel_homo_RD_20.LargerThan100.vcf > Pindel_QC_homo_RD20_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_homo_RD20_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD20_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD20_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD20_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD20_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD20_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD20_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD20_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD20_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD20_LargerThan100.INV.Mappable.TRAFree.bed

vcf.size.filter.py -i Pindel_homo_RD_30.vcf --size 100
grep -v LowQual Pindel_homo_RD_30.LargerThan100.vcf > Pindel_QC_homo_RD30_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_homo_RD30_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD30_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD30_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD30_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD30_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD30_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD30_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD30_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD30_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD30_LargerThan100.INV.Mappable.TRAFree.bed

vcf.size.filter.py -i Pindel_homo_RD_40.vcf --size 100
grep -v LowQual Pindel_homo_RD_40.LargerThan100.vcf > Pindel_QC_homo_RD40_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_homo_RD40_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD40_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD40_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD40_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD40_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD40_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD40_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD40_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD40_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD40_LargerThan100.INV.Mappable.TRAFree.bed

vcf.size.filter.py -i Pindel_homo_RD_50.vcf --size 100
grep -v LowQual Pindel_homo_RD_50.LargerThan100.vcf > Pindel_QC_homo_RD50_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_homo_RD50_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD50_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD50_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_homo_RD50_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD50_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD50_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.homo.FussyJunc/sv_rec/homo.homo.TRA.rec --input Pindel_QC_homo_RD50_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD50_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD50_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_homo_RD50_LargerThan100.INV.Mappable.TRAFree.bed

ln -s Pindel_QC_homo_RD10_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_homo_RD10_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
ln -s Pindel_QC_homo_RD20_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_homo_RD20_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
ln -s Pindel_QC_homo_RD30_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_homo_RD30_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
ln -s Pindel_QC_homo_RD40_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_homo_RD40_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
ln -s Pindel_QC_homo_RD50_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_homo_RD50_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
```
Process ref SV for comparison:
```
SV.Output.Process.py Mappable-Control --input homo.homo.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input homo.homo.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input homo.homo.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input homo.homo.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input homo.homo.DEL.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input homo.homo.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input homo.homo.DUP_TANDEM.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input homo.homo.INV.Mappable.bed

ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Delly_QC_homo_RD10_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Delly_QC_homo_RD10_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_homo_RD10_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Delly_QC_homo_RD10_INV.INV.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Delly_QC_homo_RD20_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Delly_QC_homo_RD20_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_homo_RD20_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Delly_QC_homo_RD20_INV.INV.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Delly_QC_homo_RD30_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Delly_QC_homo_RD30_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_homo_RD30_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Delly_QC_homo_RD30_INV.INV.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Delly_QC_homo_RD40_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Delly_QC_homo_RD40_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_homo_RD40_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Delly_QC_homo_RD40_INV.INV.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Delly_QC_homo_RD50_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Delly_QC_homo_RD50_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_homo_RD50_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Delly_QC_homo_RD50_INV.INV.Mappable.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD10_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD10_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD10_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD10_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD20_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD20_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD20_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD20_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD30_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD30_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD30_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD30_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD40_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD40_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD40_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD40_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD50_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD50_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD50_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD50_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD10_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD10_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD10_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD10_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD20_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD20_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD20_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD20_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD30_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD30_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD30_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD30_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD40_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD40_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD40_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD40_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD50_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD50_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD50_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed Pindel_QC_homo_RD50_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD10_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD10_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD10_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD10_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD20_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD20_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD20_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD20_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD30_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD30_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD30_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD30_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD40_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD40_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD40_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD40_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DEL.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD50_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD50_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD50_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s homo.homo.INV.Mappable.min100.max1000000000.bed SVelter_QC_homo_RD50_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
```
######Run comparison algorithm:
```
Produce.Pseudo.ROC.stats.py --path_ref rec_sv_TRAFree_rec2_20150903/ --path_in alt_sv_TRAFree_rec2_20150903/ --appdix .Mappable.TRAFree.min100.max1000000000.bed
```
>



#####To process **Simulated Simple Het** events:
<
######SVelter:
```
change.vcf.name.py -p ./

grep -v LowQual het_RD10_sorted.vcf > SVelter_QC_het_RD10_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_het_RD10_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD10_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD10_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD10_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD10_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD10_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD10_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD10_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD10_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD10_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD10_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD10_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD10_sorted.DUP.Mappable.TRAFree.bed

grep -v LowQual het_RD20_sorted.vcf > SVelter_QC_het_RD20_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_het_RD20_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD20_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD20_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD20_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD20_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD20_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD20_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD20_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD20_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD20_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD20_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD20_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD20_sorted.DUP.Mappable.TRAFree.bed

grep -v LowQual het_RD30_sorted.vcf > SVelter_QC_het_RD30_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_het_RD30_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD30_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD30_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD30_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD30_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD30_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD30_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD30_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD30_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD30_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD30_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD30_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD30_sorted.DUP.Mappable.TRAFree.bed

grep -v LowQual het_RD40_sorted.vcf > SVelter_QC_het_RD40_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_het_RD40_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD40_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD40_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD40_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD40_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD40_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD40_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD40_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD40_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD40_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD40_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD40_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD40_sorted.DUP.Mappable.TRAFree.bed

grep -v LowQual het_RD50_sorted.vcf > SVelter_QC_het_RD50_sorted.vcf
SV.Output.Process.py vcf --input SVelter_QC_het_RD50_sorted.vcf 
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD50_sorted.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD50_sorted.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD50_sorted.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input SVelter_QC_het_RD50_sorted.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD50_sorted.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD50_sorted.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD50_sorted.INV.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input SVelter_QC_het_RD50_sorted.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD50_sorted.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD50_sorted.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD50_sorted.INV.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input SVelter_QC_het_RD50_sorted.DUP.Mappable.TRAFree.bed
```
######Delly
```
change.vcf.name.py -p ./

grep -v LowQual het_RD10_DEL.vcf > Delly_QC_het_RD10_DEL.vcf
grep -v LowQual het_RD10_DUP.vcf > Delly_QC_het_RD10_DUP.vcf
grep -v LowQual het_RD10_INV.vcf > Delly_QC_het_RD10_INV.vcf
grep -v LowQual het_RD10_TRA.vcf > Delly_QC_het_RD10_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_het_RD10_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD10_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD10_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD10_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD10_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD10_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD10_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD10_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD10_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD10_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD10_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD10_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD10_INV.INV.Mappable.TRAFree.bed

grep -v LowQual het_RD20_DEL.vcf > Delly_QC_het_RD20_DEL.vcf
grep -v LowQual het_RD20_DUP.vcf > Delly_QC_het_RD20_DUP.vcf
grep -v LowQual het_RD20_INV.vcf > Delly_QC_het_RD20_INV.vcf
grep -v LowQual het_RD20_TRA.vcf > Delly_QC_het_RD20_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_het_RD20_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD20_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD20_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD20_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD20_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD20_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD20_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD20_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD20_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD20_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD20_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD20_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD20_INV.INV.Mappable.TRAFree.bed

grep -v LowQual het_RD30_DEL.vcf > Delly_QC_het_RD30_DEL.vcf
grep -v LowQual het_RD30_DUP.vcf > Delly_QC_het_RD30_DUP.vcf
grep -v LowQual het_RD30_INV.vcf > Delly_QC_het_RD30_INV.vcf
grep -v LowQual het_RD30_TRA.vcf > Delly_QC_het_RD30_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_het_RD30_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD30_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD30_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD30_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD30_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD30_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD30_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD30_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD30_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD30_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD30_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD30_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD30_INV.INV.Mappable.TRAFree.bed

grep -v LowQual het_RD40_DEL.vcf > Delly_QC_het_RD40_DEL.vcf
grep -v LowQual het_RD40_DUP.vcf > Delly_QC_het_RD40_DUP.vcf
grep -v LowQual het_RD40_INV.vcf > Delly_QC_het_RD40_INV.vcf
grep -v LowQual het_RD40_TRA.vcf > Delly_QC_het_RD40_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_het_RD40_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD40_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD40_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD40_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD40_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD40_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD40_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD40_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD40_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD40_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD40_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD40_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD40_INV.INV.Mappable.TRAFree.bed

grep -v LowQual het_RD50_DEL.vcf > Delly_QC_het_RD50_DEL.vcf
grep -v LowQual het_RD50_DUP.vcf > Delly_QC_het_RD50_DUP.vcf
grep -v LowQual het_RD50_INV.vcf > Delly_QC_het_RD50_INV.vcf
grep -v LowQual het_RD50_TRA.vcf > Delly_QC_het_RD50_TRA.vcf
SV.Output.Process.py vcf --input Delly_QC_het_RD50_DEL.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD50_DUP.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD50_INV.vcf  
SV.Output.Process.py vcf --input Delly_QC_het_RD50_TRA.vcf  
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD50_DEL.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD50_DUP.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Delly_QC_het_RD50_INV.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD50_DEL.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD50_DUP.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Delly_QC_het_RD50_INV.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD50_DEL.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD50_DUP.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Delly_QC_het_RD50_INV.INV.Mappable.TRAFree.bed

ln -s Delly_QC_het_RD10_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_het_RD10_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Delly_QC_het_RD20_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_het_RD20_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Delly_QC_het_RD30_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_het_RD30_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Delly_QC_het_RD40_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_het_RD40_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Delly_QC_het_RD50_DUP.DUP.Mappable.TRAFree.min100.max1000000000.bed Delly_QC_het_RD50_DUP_TANDEM.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
```
######Lumpy
```
change.bedpe.name.py -p ./

mv het_RD10_pesr.bedpe Lumpy_QC_RD10_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD10_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD10_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD10_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD10_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD10_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD10_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD10_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD10_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD10_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD10_pesr.INV.Mappable.TRAFree.bed

mv het_RD20_pesr.bedpe Lumpy_QC_RD20_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD20_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD20_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD20_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD20_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD20_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD20_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD20_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD20_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD20_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD20_pesr.INV.Mappable.TRAFree.bed

mv het_RD30_pesr.bedpe Lumpy_QC_RD30_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD30_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta 
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD30_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD30_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD30_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD30_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD30_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD30_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD30_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD30_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD30_pesr.INV.Mappable.TRAFree.bed

mv het_RD40_pesr.bedpe Lumpy_QC_RD40_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD40_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta 
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD40_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD40_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD40_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD40_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD40_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD40_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD40_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD40_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD40_pesr.INV.Mappable.TRAFree.bed

mv het_RD50_pesr.bedpe Lumpy_QC_RD50_pesr.bedpe
SV.Output.Process.py bedpe --input Lumpy_QC_RD50_pesr.bedpe --reference /scratch/remills_flux/xuefzhao/reference/GRCh37/human_g1k_v37.fasta 
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD50_pesr.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD50_pesr.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Lumpy_QC_RD50_pesr.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD50_pesr.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD50_pesr.DUP.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Lumpy_QC_RD50_pesr.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD50_pesr.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD50_pesr.DUP.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Lumpy_QC_RD50_pesr.INV.Mappable.TRAFree.bed

ln -s Lumpy_QC_RD10_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD10_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Lumpy_QC_RD20_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD20_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Lumpy_QC_RD30_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD30_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Lumpy_QC_RD40_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD40_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
ln -s Lumpy_QC_RD50_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed Lumpy_QC_RD50_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed
```
######Pindel
```
change.vcf.name.py -p ./

vcf.size.filter.py -i Pindel_het_RD_10.vcf --size 100
grep -v LowQual Pindel_het_RD_10.LargerThan100.vcf > Pindel_QC_het_RD10_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_het_RD10_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD10_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD10_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD10_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD10_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD10_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD10_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD10_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD10_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD10_LargerThan100.INV.Mappable.TRAFree.bed

vcf.size.filter.py -i Pindel_het_RD_20.vcf --size 100
grep -v LowQual Pindel_het_RD_20.LargerThan100.vcf > Pindel_QC_het_RD20_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_het_RD20_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD20_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD20_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD20_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD20_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD20_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD20_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD20_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD20_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD20_LargerThan100.INV.Mappable.TRAFree.bed

vcf.size.filter.py -i Pindel_het_RD_30.vcf --size 100
grep -v LowQual Pindel_het_RD_30.LargerThan100.vcf > Pindel_QC_het_RD30_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_het_RD30_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD30_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD30_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD30_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD30_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD30_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD30_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD30_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD30_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD30_LargerThan100.INV.Mappable.TRAFree.bed

vcf.size.filter.py -i Pindel_het_RD_40.vcf --size 100
grep -v LowQual Pindel_het_RD_40.LargerThan100.vcf > Pindel_QC_het_RD40_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_het_RD40_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD40_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD40_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD40_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD40_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD40_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD40_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD40_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD40_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD40_LargerThan100.INV.Mappable.TRAFree.bed

vcf.size.filter.py -i Pindel_het_RD_50.vcf --size 100
grep -v LowQual Pindel_het_RD_50.LargerThan100.vcf > Pindel_QC_het_RD50_LargerThan100.vcf
SV.Output.Process.py vcf --input Pindel_QC_het_RD50_LargerThan100.vcf
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD50_LargerThan100.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD50_LargerThan100.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input Pindel_QC_het_RD50_LargerThan100.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD50_LargerThan100.DEL.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD50_LargerThan100.DUP_TANDEM.Mappable.bed
SV.Output.Process.py TRA-Control --TRA-rec /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Simulate.het.FussyJunc/sv_rec/het.het.TRA.rec --input Pindel_QC_het_RD50_LargerThan100.INV.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD50_LargerThan100.DEL.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD50_LargerThan100.DUP_TANDEM.Mappable.TRAFree.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input Pindel_QC_het_RD50_LargerThan100.INV.Mappable.TRAFree.bed

ln -s Pindel_QC_het_RD10_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_het_RD10_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
ln -s Pindel_QC_het_RD20_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_het_RD20_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
ln -s Pindel_QC_het_RD30_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_het_RD30_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
ln -s Pindel_QC_het_RD40_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_het_RD40_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
ln -s Pindel_QC_het_RD50_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed Pindel_QC_het_RD50_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed
```


Process ref SV for comparison:
```
SV.Output.Process.py Mappable-Control --input het.het.DEL.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input het.het.DUP.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input het.het.DUP_TANDEM.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Mappable-Control --input het.het.INV.bed --ref-prefix /scratch/remills_flux/xuefzhao/reference/GRCh37_decoy/human_g1k_v37_decoy.Mappable
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input het.het.DEL.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input het.het.DUP.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input het.het.DUP_TANDEM.Mappable.bed
SV.Output.Process.py Size-Control --min-size 100 --max-size 1000000000 --input het.het.INV.Mappable.bed

ln -s het.het.DEL.Mappable.min100.max1000000000.bed Delly_QC_het_RD10_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Delly_QC_het_RD10_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_het_RD10_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Delly_QC_het_RD10_INV.INV.Mappable.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Delly_QC_het_RD20_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Delly_QC_het_RD20_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_het_RD20_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Delly_QC_het_RD20_INV.INV.Mappable.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Delly_QC_het_RD30_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Delly_QC_het_RD30_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_het_RD30_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Delly_QC_het_RD30_INV.INV.Mappable.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Delly_QC_het_RD40_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Delly_QC_het_RD40_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_het_RD40_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Delly_QC_het_RD40_INV.INV.Mappable.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Delly_QC_het_RD50_DEL.DEL.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Delly_QC_het_RD50_DUP.DUP.Mappable.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Delly_QC_het_RD50_DUP_TANDEM.DUP_TANDEM.Mappable.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Delly_QC_het_RD50_INV.INV.Mappable.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD10_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD10_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD10_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD10_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD20_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD20_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD20_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD20_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD30_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD30_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD30_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD30_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD40_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD40_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD40_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD40_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Lumpy_QC_RD50_pesr.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Lumpy_QC_RD50_pesr.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Lumpy_QC_RD50_pesr.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Lumpy_QC_RD50_pesr.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Pindel_QC_het_RD10_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Pindel_QC_het_RD10_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_het_RD10_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Pindel_QC_het_RD10_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Pindel_QC_het_RD20_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Pindel_QC_het_RD20_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_het_RD20_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Pindel_QC_het_RD20_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Pindel_QC_het_RD30_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Pindel_QC_het_RD30_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_het_RD30_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Pindel_QC_het_RD30_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Pindel_QC_het_RD40_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Pindel_QC_het_RD40_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_het_RD40_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Pindel_QC_het_RD40_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed Pindel_QC_het_RD50_LargerThan100.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed Pindel_QC_het_RD50_LargerThan100.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed Pindel_QC_het_RD50_LargerThan100.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed Pindel_QC_het_RD50_LargerThan100.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed SVelter_QC_het_RD10_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed SVelter_QC_het_RD10_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_het_RD10_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed SVelter_QC_het_RD10_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed SVelter_QC_het_RD20_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed SVelter_QC_het_RD20_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_het_RD20_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed SVelter_QC_het_RD20_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed SVelter_QC_het_RD30_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed SVelter_QC_het_RD30_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_het_RD30_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed SVelter_QC_het_RD30_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed SVelter_QC_het_RD40_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed SVelter_QC_het_RD40_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_het_RD40_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed SVelter_QC_het_RD40_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DEL.Mappable.min100.max1000000000.bed SVelter_QC_het_RD50_sorted.DEL.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP.Mappable.min100.max1000000000.bed SVelter_QC_het_RD50_sorted.DUP.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.DUP_TANDEM.Mappable.min100.max1000000000.bed SVelter_QC_het_RD50_sorted.DUP_TANDEM.Mappable.TRAFree.min100.max1000000000.bed 
ln -s het.het.INV.Mappable.min100.max1000000000.bed SVelter_QC_het_RD50_sorted.INV.Mappable.TRAFree.min100.max1000000000.bed 
```
######Run comparison algorithm:
```
Produce.Pseudo.ROC.stats.py --path_ref rec_sv_TRAFree_rec2_20150903/ --path_in alt_sv_TRAFree_rec2_20150903/ --appdix .Mappable.TRAFree.min100.max1000000000.bed
```
>
