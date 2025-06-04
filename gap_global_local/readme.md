Run check_networks(L, nb_modules) to check EM networks with L nodes and nb_modules unknown modules for a potential EM network which would be generically locally identifiable, but not generically globally identifiable.

This checks all EM networks with L nodes, nb_modules unknown modules and no known modules, for which nB * nC = nb_modules (where nB is the number of excitations, and nC the number of measurements), except:
1. The ones with full excitation or full measurement
2. The disconnected ones
3. The ones with self-loops

Description of each function:

- **check_networks.m**
is the main function

- **solve_CTB_CTbisB.m**
is an auxiliary function that computes the equation CTB = CT'B symbolically, which yields the generic global identifiability of the network

- **is_connected.m**
is an auxiliary function that checks whether the graph is connected
