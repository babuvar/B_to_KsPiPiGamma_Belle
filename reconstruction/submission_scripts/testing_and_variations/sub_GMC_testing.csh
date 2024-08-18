#! /bin/tcsh -f

foreach Stream (0)

foreach Type (charged mixed charm uds)

#                               run_st  r_end   exp     type     mc-stream 
bsub -q s ./script_GMC_4s.csh 	64	115	19	$Type 	 $Stream

# break to indicate change in Stream-argument-format from this experiment onward

bsub -q s ./script_GMC_4s.csh 	589	620	37	$Type 	 $Stream

end

end


