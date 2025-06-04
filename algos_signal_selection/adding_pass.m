% Randomly adds signals until rank(K) = nbunknowns

function [excited_idx, measured_idx] = adding_pass(graph, unknowns, ...
    Kfull, excited_idx_0, measured_idx_0)

L = length(graph);
nbunknowns = sum(unknowns, 'all');

% Initialization
excited_idx = excited_idx_0;
excited = zeros(L,1);
excited(excited_idx) = 1;
measured_idx = measured_idx_0;
measured = zeros(L,1);
measured(measured_idx) = 1;
nB = length(excited_idx);
nC = length(measured_idx);
indices = indexing_K(L, excited_idx, measured_idx);

while rank(Kfull(indices,:)) < nbunknowns
    nb_not_excited = nnz(~ excited);
    nb_not_measured = nnz(~ measured);
    nb_not = nb_not_excited + nb_not_measured;
    if rand * nb_not < nb_not_excited
        % Add excitation
        not_excited_idx = find(~ excited);
        new_excitation = ceil(rand * nb_not_excited);
        new_excitation_idx = not_excited_idx(new_excitation);
        excited(new_excitation_idx) = 1;
        excited_idx = [excited_idx; new_excitation_idx];
        nB = nB + 1;
    else
        % Add measurement
        not_measured_idx = find(~ measured);
        new_measurement = ceil(rand * nb_not_measured);
        new_measurement_idx = not_measured_idx(new_measurement);
        measured(new_measurement_idx) = 1;
        measured_idx = [measured_idx; new_measurement_idx];
        nC = nC + 1;
    end
    indices = indexing_K(L, excited_idx, measured_idx);
end