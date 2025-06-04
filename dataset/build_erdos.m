% Builds an Erdős-Rényi graph of L nodes and 3*L/2 edges, without self-loops.
% Makes 10 selections of L/2 unknown modules among those edges.
% Stores the networks in a .m file
% If display == 1, diplays each network generated.
% The random seed is fixed, so each call to this function will output the
% same networks.

function build_erdos(L, display)
close all;

ratio = 1/3;
ntrials = 10;

% Parameters
density = 1.5;
nb_edges = round(density * L);

rng(1);
graph_idx_vect_nodiag = randperm(L*(L-1), nb_edges);
graph_vect_nodiag = zeros(L*(L-1), 1);
graph_vect_nodiag(graph_idx_vect_nodiag) = 1;
graph = zeros(L,L);
count = 1;
for i = 1:L
    for j = 1:L
        % Avoiding self-loops
        if i ~= j
            graph(i,j) = graph_vect_nodiag(count);
            count = count + 1;
        end
    end
end


% Generating unknown data for ntrials samples
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

topology = "erdos";
filename = topology + string(L) + ".mat";

graph = sparse(graph);
%unknowns_tensor = sparse(unknowns_tensor); Matlab doesnt support sparse 3D
save(filename, 'graph', 'unknowns_tensor', 'topology', 'L');
excited = zeros(L,1);
measured = zeros (L,1);

% Display networks with unknowns
if display
    for i = 1:ntrials
        main_identifiable("Erdős-Rényi graph", graph, unknowns_tensor(:,:,i), excited, measured, 1, 0);
    end
end

end