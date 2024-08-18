#! /bin/tcsh -f



foreach i (b_kspipigam_Sm0p6_1M  b_kspipigam_Sm0p4_1M   b_kspipigam_Sm0p2_1M  b_kspipigam_S0p0_1M  b_kspipigam_S0p2_1M  b_kspipigam_S0p4_1M  b_kspipigam_S0p6_1M)

foreach j (`seq 1 38`)

#bsub -q l ./script_SMC.csh $i $j
bsub -q s ./script_SMC.csh $i $j



end

end




