% Builds a Lattice graph of L nodes
% Makes 10 selections of 1/3 of the edges as unknown modules.
% Stores the networks in a .m file
% If display == 1, diplays each network generated.
% The random seed is fixed, so each call to this function will output the
% same networks.

function build_lattice(rows, cols, display)
close all;

% Total number of nodes
L = rows * cols;
ratio = 1/3;
ntrials = 10;

% Create a sparse adjacency matrix
graph = zeros(L, L);

% Construct the lattice adjacency matrix
seed = 1;
for r = 1:rows
    for c = 1:cols
        node = (r-1) * cols + c;  % Current node index

        % Connect to the right neighbor
        if c < cols
            rightNeighbor = node + 1;
            rng(seed);
            if rand < 0.5
                graph(node, rightNeighbor) = 1;
            else
                graph(rightNeighbor, node) = 1;
            end
            seed = seed + 1;
        end

        % Connect to the bottom neighbor
        if r < rows
            bottomNeighbor = node + cols;
            rng(seed);
            if rand < 0.5
                graph(node, bottomNeighbor) = 1;
            else
                graph(bottomNeighbor, node) = 1;
            end
            seed = seed + 1;
        end
    end
end

% Generating unknown data for ntrials samples
nb_edges = nnz(graph);
m = round(nb_edges * ratio); % nb unknowns
[row_edges, col_edges] = find(graph);
unknowns_tensor = zeros(L,L,ntrials);
for i = 1:ntrials
    rng(10*i); % Fix seed for each i
    % Randomly selects m unknowns among the edges
    random_idx = randperm(nb_edges, m);
    for j = 1:m
        unknowns_tensor(row_edges(random_idx(j)), col_edges(random_idx(j)), i) = 1;
    end
end

topology = "lattice";
filename = topology + string(rows) + "x" + string(cols) + ".mat";

graph = sparse(graph);
%unknowns_tensor = sparse(unknowns_tensor); Matlab doesnt support sparse 3D
save(filename, 'graph', 'unknowns_tensor', 'topology', 'L');
excited = zeros(L,1);
measured = zeros (L,1);

% Display networks with unknowns
if display
    for i = 1:ntrials
        main_identifiable("Lattice graph", graph, unknowns_tensor(:,:,i), excited, measured, 1, 0);
    end
end

end