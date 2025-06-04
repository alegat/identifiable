% From initial sets of excited and measured nodes, tests if generic local
% identifiability is conserved under each excitation and measurement removal,
% 1 by 1.
% This is done in a uniform random order npass times, and the best output
% is returned.

function [excited_idx_best, measured_idx_best] = pruning_multipass(graph, unknowns, ...
    Kfull, excited_idx_0, measured_idx_0, npass)

L = length(graph);
nbunknowns = sum(unknowns, 'all');
nbsignals_bound = compute_bound(nbunknowns);

excited_idx_best = excited_idx_0;
measured_idx_best = measured_idx_0;
nbsignals_0 = length(excited_idx_0) + length(measured_idx_0);
nbsignals_best = nbsignals_0;

if nbsignals_0 > nbsignals_bound

    for pass = 1:npass

        % Pruning pass
        list = 1:nbsignals_0;
        shuffled_list = list(randperm(length(list)));
        nb_exc_0 = length(excited_idx_0);
        excited_idx = excited_idx_0;
        measured_idx = measured_idx_0;

        for i = 1:nbsignals_0
            node_idx = shuffled_list(i);
            excited_idx_new = excited_idx;
            measured_idx_new = measured_idx;
            if node_idx <= nb_exc_0
                excited_idx_new(excited_idx_new == excited_idx_0(node_idx)) = [];
            else
                measured_idx_new(measured_idx_new == measured_idx_0(node_idx - nb_exc_0)) = [];
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

end