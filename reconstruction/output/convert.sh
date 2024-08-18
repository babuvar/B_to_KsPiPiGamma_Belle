#!/usr/bin/sh

# Usage : ./convert.sh SMC
case $1 in

SMC)
cd hbook_SMC
for file in *.hbk
do
h2root $file
done
mv *.root ../root_SMC
;;

# Usage : ./convert.sh GMC 2(stream no.)
GMC)
cd hbook_GMC$2
for file in *.hbk
do
h2root $file
done
mv *.root ../root_GMC$2
;;

# Usage : ./convert.sh data
data)
cd hbook_data
for file in *.hbk
do
h2root $file
done
mv *.root ../root_data
;;

# Usage : ./convert.sh rare
rare)
cd hbook_MCrare
for file in *.hbk
do
h2root $file
done
mv *.root ../root_MCrare
;;


esac

