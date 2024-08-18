#! /bin/tcsh -f

setenv subdir s0

foreach exp ( e07  e09  e11  e13  e15  e17  e19a  e19b  e21  e23  e25  e27  e61  e63  e65 )

bsub -q s ./script_MCrare.csh	$subdir	     $exp

end





setenv subdir s1

foreach exp (e31  e33  e35  e37a  e37b  e39  e41a  e41b  e43  e45a  e45b  e47  e49  e51  e55 )

bsub -q s ./script_MCrare.csh   $subdir      $exp

end


