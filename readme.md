Those 4 functions compute, print and plot the generic identifiability of 6 small networks taken from [HGB18] and [CSV19] (see reference below).

# test_identifiable.m
is the only file you need to play with. We define the graphs to test through their adjacency matrix, the excited and measured nodes. We then call main_identifiable.

## main_identifiable.m
takes data as argument, passes to identifiable.m that returns the solution, prints it and calls plot_identifiable.m for a visualisation.

### identifiable.m 
does all the computation. Note that you can choose the number of random samples (although 1 should be enough), or opt for symbolic computation (slower).

### plot_identifiable.m
plots the solution.

[HGB18] Hendrickx, Julien M., Michel Gevers, and Alexandre S. Bazanella. "Identifiability of dynamical networks with partial node measurements." IEEE Transactions on Automatic Control 64.6 (2018): 2240-2253.

[CSV19] Cheng, Xiaodong, Shengling Shi, and Paul MJ Van den Hof. "Allocation of excitation signals for generic identifiability of linear dynamic networks." arXiv preprint arXiv:1910.04525v1 (2019).
