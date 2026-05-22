#!/bin/bash

# Initialize CSV files
for opt_level in O0 O1 O2 O3 O4; do
    echo "nodes,density,exp,delta,compressionRatio,executionTime" > "cpgr_results_${opt_level}.csv"
    echo "nodes,density,exp,delta,compressionRatio,executionTime" > "fm_results_${opt_level}.csv"

    echo "nodes,total_nodes,density,exp,delta,maximumFlow,run_time,total_run_time" > "bipartite_dinics_results_${opt_level}.csv"
    echo "nodes,total_nodes,density,exp,delta,maximumFlow,run_time,total_run_time" > "tripartite_dinics_results_${opt_level}.csv"
done

algorithm="cpgr"
for node in 32 64 128 256 512 1024 2048 4096 8192 16384 32768; do
    for density in 80 85 90 95 98; do
        for exp in {1..10}; do
            for delta in 0.5 0.6 0.7 0.8 0.9 1; do
                # Run cpgr algorithms
                for opt_level in O0 O1 O2 O3 O4; do
                    ./cpgr_$opt_level "$node" "$density" "$exp" "$delta" >> "cpgr_results_${opt_level}.csv"
                done

                # Conditionally run fm algorithm for small graphs
                if [ "$node" -le 128 ]; then
                    for opt_level in O0 O1 O2 O3 O4; do
                        ./fm_$opt_level "$node" "$density" "$exp" "$delta" >> "fm_results_${opt_level}.csv"
                    done
                fi
                
                # executing the dinitz algorithm with different optimization level for given and compressed graph
                for opt_level in o0 o1 o2 o3 o4; do
                    ./dinics_bi_$opt_level "$node" "$density" "$exp" "$delta" >> "bipartite_dinics_results_${opt_level}.csv"
                    ./dinics_tri_$opt_level "$node" "$density" "$exp" "$delta" "$algorithm" >> "tripartite_dinics_results_${opt_level}.csv"
                done

                echo "$node $density $exp"
            done
        done
    done
done
