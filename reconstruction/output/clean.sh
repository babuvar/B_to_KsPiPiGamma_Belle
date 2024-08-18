#!/usr/bin/sh

# Usage : ./clean.sh SMC root/hbook 3(0...5)
case $1 in

SMC)
case $2 in
root)
cd root_SMC
rm -f *.root
;;

hbook)
cd hbook_SMC
rm -f *.hbk
;;
esac
;;

# Usage : ./clean.sh GMC root/hbook 3(0...5)
GMC)
case $2 in
root)
cd root_GMC$3
rm -f *.root
;;

hbook)
cd hbook_GMC$3
rm -f *.hbk
;;
esac
;;

# Usage : ./clean.sh data root/hbook
data)
case $2 in
root)
cd root_data
rm -f *.root
;;

hbook)
cd hbook_data
rm -f *.hbk
;;
esac
;;

# Usage : ./clean.sh rare root/hbook 
rare)
case $2 in
root)
cd root_MCrare
rm -f *.root
;;

hbook)
cd hbook_MCrare
rm -f *.hbk
;;
esac
;;

esac

