


#for i in {1..500};
for i in {1..2500};
do
    j=$(expr $i - 1 )
    bsub -q s ./run_toyFit.sh $j
done

