cd /gpfs/home/belle/varghese/IPHC/B_Kspipigamma_belle/tatami_fits/alternate_way

#for i in {1..100};
#for i in {1..1000};
#for i in {1..5};
#do
i=1

for smc_sample in S0p2 S0p4 S0p6 Sm0p2 Sm0p4 Sm0p6 S0p0
#for smc_sample in S0p2

do
    j=$(expr $i - 1 )
    #bsub -q s ./run_linearityFit_v2.sh $j $smc_sample
    #bsub -q s ./run_linearityFit_v2_3DShapeFix.sh $j $smc_sample
    #bsub -q s ./run_linearityFit_v2.sh $j $smc_sample
    #bsub -q s python3 run_linearityFit_v2.py $j $smc_sample
    python3 run_linearityFit_v2.py $j $smc_sample

done

#done
