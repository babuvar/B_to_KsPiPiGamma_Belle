export nToys_per_job=10

#loop over different SMC parent-samples
#for SMC_sample in S0p2 S0p4 S0p6 Sm0p0 Sm0p2 Sm0p4 Sm0p6 
#for SMC_sample in S0p6 Sm0p6
for SMC_sample in Sm0p6 Sm0p4 Sm0p2 S0p0 S0p2 S0p4 S0p6

do

#Submit 100 jobs for each SMC sample
for i in {1..100};
do
job=$(expr $i - 1 )

bsub -q s python3 make_bootstrapped_datasets_smalljob_onlySignal.py $nToys_per_job $job $SMC_sample

done

done

