#!/bin/bash

echo "nodes,density,exp,delta,compression_ratio,execution_time"  > fm_results.csv
for node in  32 64 128 
do
    	for density in 80 85 90 95 98 
    	do
        	for exp in 1 2 3 4 5 6 7 8 9 10
        	do
	 		for delta in 0.5 0.6 0.7 0.8 0.9 1
    			do
				./fm  $node $density $exp $delta >> fm_results.csv
            			echo $node $density $exp
	       		done
        	done
    	done
done

