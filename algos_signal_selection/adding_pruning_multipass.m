% From initial sets of excited and measured nodes, randomly adds
% excitations and measurements if not generically locally identifiable.
% Once identifiability has been reached, it tests if generic local
% identifiability is conserved under each excitation and measurement removal,
% 1 by 1.
% This is done in a uniform random order npass times, and the best output
% is returned.

function [excited_idx_best, measured_idx_best] = adding_pruning_multipass...
    (graph, unknowns, Kfull, excited_idx_0, measured_idx_0, npass)

L = length(graph);
nbunknowns = sum(unknowns, 'all');
nbsignals_bound = compute_bound(nbunknowns);

% Initialization
excited_idx_best = ones(L,1);
measured_idx_best = ones(L,1);
nbsignals_best = 2*L;

for pass = 1:npass
    excited_idx = excited_idx_0;
    measured_idx = measured_idx_0;
    indices = indexing_K(L, excited_idx, measured_idx);

    if rank(Kfull(indices, :)) < nbunknowns
        % Adding pass
        [excited_idx, measured_idx] = adding_pass(graph, unknowns, ...
            Kfull, excited_idx, measured_idx);
    end

    excited_idx_1 = excited_idx;
    measured_idx_1 = measured_idx;
    nb_exc_1 = length(excited_idx_1);
    nbsignals_1 = nb_exc_1 + length(measured_idx_1);
    list = 1:nbsignals_1;
    shuffled_list = list(randperm(length(list)));

    % Pruning pass
    for i = list
        node_idx = shuffled_list(i);
        excited_idx_new = excited_idx;
        measured_idx_new = measured_idx;
        if node_idx <= nb_exc_1
            excited_idx_new(excited_idx_new == excited_idx_1(node_idx)) = [];
        else
            measured_idx_new(measured_idx_new == measured_idx_1(node_idx - nb_exc_1)) = [];
        end

        indices = indexing_K(L, excited_idx_new, measured_idx_new);
        % Generic local identifiability test
        if rank(Kfull(indices,:)) == nbunknowns
            excited_idx = excited_idx_new;
            measured_idx = measured_idx_new;
        end
    end

    % Optimality test
    nbsignals = length(excited_idx) + length(measured_idx);
    if nbsignals <= nbsignals_bound
        % Cannot do better
        excited_idx_best = excited_idx;
        measured_idx_best = measured_idx;
        break;
    elseif nbsignals < nbsignals_best
        nbsignals_best = nbsignals;
        excited_idx_best = excited_idx;
        measured_idx_best = measured_idx;
    end
end

end