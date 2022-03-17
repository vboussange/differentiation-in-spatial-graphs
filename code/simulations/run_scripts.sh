#########################################################################################
# This scripts launches in parallel 5 replicates of the same script with different seeds
# It needs to be placed in the same folder as the script
# The environment should be at one level above
#########################################################################################
#!/bin/bash
namesim=$1 #'neutralmarker_discrete_D2_D1_varying_mu_1e0_a_1e0'
mkdir $namesim

for i in {1..5}
do
	echo "lauching script $namesim for seed $i"
	($HOME/utils/julia-1.6.1/bin/julia --project=. --threads 6 "$namesim.jl" $i &> "${namesim}/seed${i}.out") &
done
wait
echo "computation over"
date