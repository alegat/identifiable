% Builds a Watts–Strogatz graph of L nodes
% Makes 10 selections of 1/3 of the edges as unknown modules.
% Stores the networks in a .m file
% If display == 1, diplays each network generated.
% The random seed is fixed, so each call to this function will output the
% same networks.

function build_watts(L, display)
close all;

ratio = 1/3;
ntrials = 10;
seed = 1;

% Parameters
k = 4;               % Each node connects to k nearest neighbors (even number)
rewiringProb = 0.1;  % Probability of rewiring each edge

% Create a ring lattice where each node is connected to k/2 neighbors on each side
graph = zeros(L,L);    % Initialize adjacency matrix
halfK = k / 2;         % Number of neighbors on each side

% Build the regular ring lattice
for i = 1:L
    for j = 1:halfK
        neighbor = mod(i + j - 1, L) + 1;
        rng(seed);

        % Regular wiring
        if rand > rewiringProb
            seed = seed + 1;
            rng(seed);
            if rand < 0.5
                graph(i, neighbor) = 1;
            else
                graph(neighbor, i) = 1;
            end

        % Alternative wiring
        else
            % Find a new node to connect to, avoiding self-loops and duplicate edges
            newNeighbor = i;
            while newNeighbor == i || graph(i, newNeighbor) == 1
                seed = seed + 1;
                rng(seed);
                newNeighbor = randi(L);
            end
            % Add the new edge
            seed = seed + 1;
            rng(seed);
            if rand < 0.5
                graph(i, newNeighbor) = 1;
            else
                graph(newNeighbor, i) = 1;
            end
        end
        seed = seed + 1;
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

topology = "watts";
filename = topology + string(L) + ".mat";

graph = sparse(graph);
%unknowns_tensor = sparse(unknowns_tensor); Matlab doesnt support sparse 3D
save(filename, 'graph', 'unknowns_tensor', 'topology', 'L');
excited = zeros(L,1);
measured = zeros (L,1);

% Display networks with unknowns
if display
    for i = 1:ntrials
        main_identifiable("Watts–Strogatz graph", graph, unknowns_tensor(:,:,i), excited, measured, 1, 0);
    end
end