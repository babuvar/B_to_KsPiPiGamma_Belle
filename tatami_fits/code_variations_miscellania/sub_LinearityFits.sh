export Acp_A=0.0

#for i in {1..1000};
for i in {1..500};
do

#for Acp_S in -0.4 -0.2 0.0 0.2 0.4 
for Acp_S in 0.4
do
    j=$(expr $i - 1 )
    bsub -q s ./run_linearityFit.sh $j $Acp_A  $Acp_S
done

done
