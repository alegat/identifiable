% Computation of ancestors_weighed:
% ancestors_weighted(i,j) = sum_k ((estimated number of acyclic walks j->i of length k) / 2^k)
% We consider L/4 generations, i.e. k goes from 1 to L/4.
%
% Computing the exact number of acyclic walks j->i of length k is #P-complete.
% It is tempting to try to compute them with powers of the adjacency
% matrix, but if there are cycles in the graph, the result is inaccurate.
%
% This implementation is based on powers of the adjacency matrix, from
% which most cycles are removed (but some cycles are still counted).
% Experience show that even if some cycles are still counted, the results
% of the merging algorithm is satisfactory.
% This implementation is a trade-off between accuracy and efficiency.

function ancestors_weighted = compute_ancestors_weighted(graph)

L = length(graph);
depth_ancestors = round(L/4);

acyclic_walks = zeros(L, L, depth_ancestors);
acyclic_walks(:,:,1) = graph - diag(diag(graph));
ancestors_weighted = acyclic_walks(:,:,1);

for i = 2:depth_ancestors
    walks = acyclic_walks(:,:,1) * acyclic_walks(:,:,i-1);
    for j = 2:i-1
        walks = min(walks, acyclic_walks(:,:,j) * acyclic_walks(:,:,i-j));
    end
    acyclic_walks(:,:,i) = walks - diag(diag(walks));
    ancestors_weighted = ancestors_weighted + acyclic_walks(:,:,i) / 2^i;
end

end