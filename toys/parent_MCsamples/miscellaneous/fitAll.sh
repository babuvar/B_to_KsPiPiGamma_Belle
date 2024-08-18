mkdir results plots

for i in *orrect*.root
do

root -q -b -l "genRooBCPFit.C(\"$i\")"

done
