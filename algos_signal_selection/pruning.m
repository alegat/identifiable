% The Pruning Search algorithm

function [excited, measured] = pruning(graph, unknowns, npass)

addpath('../algebraic_identifiability_test')
L = length(graph);

% Random sampling
G = rand(L,L) .* graph;
T_rand = inv(eye(L) - G);
Kfull = kron_submatrix(T_rand', T_rand, unknowns);

% Start iterate
excited_idx_0 = 1:L;
measured_idx_0 = 1:L;

% Pruning passes
[excited_idx, measured_idx] = pruning_multipass(graph, unknowns, ...
    Kfull, excited_idx_0, measured_idx_0, npass);

excited = zeros(L,1); excited(excited_idx) = 1;
measured = zeros(L,1); measured(measured_idx) = 1;

end
