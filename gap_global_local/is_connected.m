% Outputs 1 if the graph is connected, 0 otherwise

function connectedness = is_connected(graph)

threshold = 1e-9;

A_nonbin = graph + graph';
A = double(A_nonbin); % undirected version of the adjacency matrix

D = diag(sum(A, 2)); % degree matrix
L = D - A; % graph laplacian

eigenvalues = sort(eig(L));  % Compute and sort eigenvalues
algebraic_connectivity = eigenvalues(2);  % Second smallest eigenvalue

connectedness = algebraic_connectivity > threshold;

end