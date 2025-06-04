% The Exhaustive Search algorithm

function [excited, measured] = exhaustive(graph, unknowns)

addpath('../algebraic_identifiability_test')
L = length(graph);
Id = eye(L);
nbunknowns = sum(unknowns(:));

% Random sampling
G = rand(L,L);
G = G .* graph;
Kfull = kron_submatrix((inv(Id-G))', inv(Id-G), unknowns);

start_it = compute_bound(nbunknowns);
for nbsignals = start_it:2*L
    nB_list = 1:(nbsignals-1);
    nC_list = nbsignals - nB_list;
    for i = 1:(nbsignals-1)
        if nB_list(i) * nC_list(i) >= nbunknowns
            nB = nB_list(i);
            nC = nC_list(i);
            combis_B = nchoosek(1:L,nB);
            combis_C = nchoosek(1:L,nC);
            [nb_combi_B, ~] = size(combis_B);
            [nb_combi_C, ~] = size(combis_C);
            for j = 1:nb_combi_B
                excited_idx = combis_B(j,:);
                for k = 1:nb_combi_C
                    measured_idx = combis_C(k,:);
                    indices = indexing_K(L, excited_idx, measured_idx);

                    % Identifiability test with rank(K)
                    if rank(Kfull(indices,:)) == nbunknowns
                        excited = zeros(L,1); excited(excited_idx) = 1;
                        measured = zeros(L,1); measured(measured_idx) = 1;
                        return
                    end
                end
            end
        end
    end
end

end
