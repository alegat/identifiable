Those four functions compute, print and plot the generic local identifiability of a few small networks.

- **test_identifiable.m**
is the only file you need to play with. We define the graphs to test through their adjacency matrix, their unknown edges, the excited and measured nodes. We then call main_identifiable.
Below is an example of input you can write in test_identifiable.m :

```
graph = [0 0 0;      % the adjacency matrix
         1 0 1;
         1 0 0];
unknowns = [0 0 0;   % the unknown edges
            1 0 0;
            1 0 0];
excited = [1 0 0]';  % the excited nodes
measured = [0 1 1]'; % the measured nodes
name = 'Example A';
nsamples = 100;
decoupled = 1;	     % Take decoupled = 1 for generic decoupled identifiability.
                     % Otherwise, it is generic local identifiability.
main_identifiable(name, graph, unknowns, excited, measured, nsamples, decoupled);
```

- **main_identifiable.m**
takes data as argument, passes to identifiable.m that returns the solution, prints it and calls plot_identifiable.m for a vizualisation.

- **identifiable.m**
does all the computation. Note that you can choose the number of random samples, or opt for symbolic computation (slower).

- **plot_identifiable.m**
provides a graphical vizualisation of the output.


