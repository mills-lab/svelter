##we also simulated simple and complex germline /  somatic SVs for comparison, as SMuFin simulations do not include enough complex events that are recently revealed popular.

For the simulation, we implemented ~15,000 simple events in tumor reference genome, out of which 2/3 were heterozygous and 1/3 were homozygous. 80% SVs in tumor reference were also implemented in germline refernece genome, leaving the rest 20% mimicing somatic events.
We also implemented 7721 complex events in tumor reference genome, of which 6177 were also impelented in germline reference.

#### Commands to produce simulated references with simple events:
```
Produce.Simulated.FussyJuncs.py case-control-simple --reference /mnt/EXT/Mills-scratch2/reference/hg19/hg19.fa --input-sim /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/simple_case_control.sim --output-prefix /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/case_control
```

#### Commands to produce simulated references with complex events:
```
Produce.Simulated.FussyJuncs.py case-control-complex --reference /mnt/EXT/Mills-scratch2/reference/hg19/hg19.fa --input-rec /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/comp_case_control.rec --input-sim /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/comp_case_control.sim --output-prefix /mnt/EXT/Mills-scratch2/Xuefang/Simulate.FussyJunc/Simulate.rerun.test.20150901/case_control
```

please files used for the simulation are also shared in this folder.


