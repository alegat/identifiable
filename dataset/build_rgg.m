% Builds a Random Geometric Graph of L nodes and 3*L/2 edges
% Makes 10 selections of L/2 unknown modules among those edges.
% Stores the networks in a .m file
% If display == 1, diplays each network generated.
% The random seed is fixed, so each call to this function will output the
% same networks.

function build_rgg(L, display)
close all;

ratio = 1/3;
ntrials = 10;
nb_edges_wanted = round(3*L/2);

% Generate random (x, y) coordinates for each node
rng(1);
x = rand(L, 1);
rng(2);
y = rand(L, 1);

distances = zeros(L,L);
graph = zeros(L,L);
% Compute pairwise distances manually
for i = 1:L
    for j = i+1:L
        distances(i,j) = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2);
    end
end

distances_nz = nonzeros(distances);
sorted_distances = sort(distances_nz(:));
radius = sorted_distances(nb_edges_wanted);

seed = 3;
for i = 1:L
    for j = i+1:L
        if distances(i,j) <= radius
            rng(seed);
            if rand < 0.5
                graph(i, j) = 1;
            else
                graph(j, i) = 1;  % Since the graph is directed
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

topology = "rgg";
filename = topology + string(L) + ".mat";

graph = sparse(graph);
%unknowns_tensor = sparse(unknowns_tensor); Matlab doesnt support sparse 3D
save(filename, 'graph', 'unknowns_tensor', 'topology', 'L');
excited = zeros(L,1);
measured = zeros (L,1);

% Display networks with unknowns
if display
    for i = 1:ntrials
        main_identifiable("RGG", graph, unknowns_tensor(:,:,i), excited, measured, 1, 0);
    end
end

end