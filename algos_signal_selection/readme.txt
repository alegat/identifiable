
Run compare_algos_performance to compare the performances of our algorithms on our dataset.

- **compare_algos_performance.m**
Runs our algorithms on all networks of the dataset with up to 50 nodes, for 10 sets of unknowns, 10 times each. Returns and saves the average performances of each algo for each network.
Make sure to run build_all.m in folder dataset before running this function, in order to generate the dataset.

- **compare_algos_performance_100.m**
Same for the networks of 100 nodes.

## Algorithms

- **exhaustive.m**
The Exhaustive Search algorithm

- **pruning.m**
The Pruning Search algorithm.
Relies on:
	- pruning_multipass.m

- **pseudotree.m**
The Pseudotree Search algorithm
Relies on:
	- adding_pruning_multipass.m
	- cheng.m
	- merge_covering.m
	- odot.m
	- odot_vec.m

- **simug.m**
The SIMUG Search algorithm
Relies on:
	- adding_pruning_multipass.m
	- dreef.m
	- merge_covering_simug.m
	- odot.m
	- otimes.m
	- odot_vec.m
	- otimes_vec.m

- **merging.m**
The Merging Search algorithm
Relies on:
	- merging_ranking_B.m
	- merging_ranking_C.m
	- merging_uniform_B.m
	- merging_uniform_C.m
	- compute_ancestors_weighted.m
	- parameters_merging.m
	- pruning_multipass.m

## Auxiliary functions

- **compute_bound.m**
returns the minimal bound on the number of signals

- **indexing_K.m**
returns the row indices for selecting a submatrix of Kfull corresponding to the sets of excitations and measurements given as arguments.
