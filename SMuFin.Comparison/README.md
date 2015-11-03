#we compared SVelter, Delly, Lumpy, Pindel and SMuFin based on the somatic events SMuFin have implemented in their paper ()

we provided here the scripts we used to to apply each algorithm and interprete the results:

####Apply Lumpy :
lumpy -mw 4 -tt 0.0 -x /scratch/remills_flux/xuefzhao/svelter/Support/Exclude.hg19.bed -pe bam_file:/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/alignment/chr22_insilico_Tumor.sorted.bam,histo_file:/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.histo,mean:499.5817,stdev:21.3189241077,read_length:80,min_non_overlap:80,discordant_z:4,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 -sr bam_file:/scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor.sorted.um.sorted.bam,back_distance:20,weight:1,id:bwa,min_mapping_threshold:20 > /scratch/remills_flux/xuefzhao/Simulation.Xuefang/Smufin/dataset/bam_file/Lumpy/chr22_insilico_Tumor
