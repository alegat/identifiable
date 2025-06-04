% The Merging Search algorithm

function [excited_best, measured_best] = merging(graph, unknowns, nrun, npass, stop)

L = length(graph);
nbunknowns = sum(unknowns(:));
[targets_unknowns, sources_unknowns] = find(unknowns);
nbsignals_bound = compute_bound(nbunknowns);
ancestors_weighted = compute_ancestors_weighted(graph);

% Initialization
excited_0 = any(unknowns,1)';
excited_idx_0 = find(excited_0);
measured_0 = any(unknowns,2);
measured_idx_0 = find(measured_0);
[sigma_C_0, sigma_B_0] = find(unknowns);
tol = 1e-9;
tried_B = zeros(L,1);
tried_C = zeros(L,1);
B_done = 0;
C_done = 0;

excited_best = excited_0;
measured_best = measured_0;
nbsignals_best = 2 * nbunknowns;

for run = 1:nrun
    % Random sampling
    G = rand(L,L) .* graph;
    T_rand = inv(eye(L) - G);
    T = abs(T_rand) > tol;
    ancestors = T - diag(diag(T));
    Kfull = kron_submatrix(T_rand', T_rand, unknowns);

    % Merging search
    switch run
        case 1
            % Ranking B version
            [excited_idx, measured_idx] = merging_ranking_B(graph, unknowns, ...
            targets_unknowns, sources_unknowns, ancestors_weighted, T_rand, Kfull, ...
            excited_idx_0, measured_idx_0, sigma_B_0, sigma_C_0, ...
            tried_B, tried_C, stop, C_done);
        case 2
            % Ranking C version 
            [excited_idx, measured_idx] = merging_ranking_C(graph, unknowns, ...
            targets_unknowns, sources_unknowns, ancestors_weighted, T_rand, Kfull, ...
            excited_idx_0, measured_idx_0, sigma_B_0, sigma_C_0, ...
            tried_B, tried_C, stop, B_done);
        otherwise
            % Uniform version
            % Toss a coin: 50% chance to start with merging_B, 50% with merging_C
            if rand < 0.5
                [excited_idx, measured_idx] = merging_uniform_B(graph, unknowns, ...
                    targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
                    excited_idx_0, measured_idx_0, sigma_B_0, sigma_C_0, stop, C_done);
            else
                [excited_idx, measured_idx] = merging_uniform_C(graph, unknowns, ...
                    targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
                    excited_idx_0, measured_idx_0, sigma_B_0, sigma_C_0, stop, B_done);
            end
    end

    % Pruning passes
    [excited_idx, measured_idx] = pruning_multipass(graph, unknowns, ...
        Kfull, excited_idx, measured_idx, npass);

    % Optimality check
    nbsignals = length(excited_idx) + length(measured_idx);
    if nbsignals <= nbsignals_bound
        % Cannot do better
        excited_best = zeros(L,1); excited_best(excited_idx) = 1;
        measured_best = zeros(L,1); measured_best(measured_idx) = 1;
        break;
    elseif nbsignals < nbsignals_best
        % Best run so far
        nbsignals_best = nbsignals;
        excited_best = zeros(L,1); excited_best(excited_idx) = 1;
        measured_best = zeros(L,1); measured_best(measured_idx) = 1;
    end
end

end

