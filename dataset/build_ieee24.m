% Encodes a Standard IEEE Power Flow Test Cases
% cfr https://icseg.iti.illinois.edu/power-cases/
% Makes 10 selections of 1/3 of the edges as unknown modules.
% Stores the networks in a .m file
% If display == 1, diplays each network generated.
% The random seed is fixed, so each call to this function will output the
% same networks.

function build_ieee24(display)
close all;

ratio = 1/3;
ntrials = 10;
L = 24;
graph = zeros(L,L);

% Branch data
graph(2,1) = 1;
graph(1,3) = 1;
graph(5,1) = 1;
graph(4,2) = 1;
graph(6,2) = 1;
graph(9,3) = 1;
graph(9,4) = 1;
graph(6,10) = 1;
graph(8,7) = 1;
graph(9,8) = 1;
graph(10,8) = 1;
graph(13,11) = 1;
graph(14,11) = 1;
graph(13,12) = 1;
graph(23,12) = 1;
graph(13,23) = 1;
graph(16,14) = 1;
graph(15,16) = 1;
graph(21,15) = 1;
graph(15,21) = 1;
graph(24,15) = 1;
graph(17,16) = 1;
graph(19,16) = 1;
graph(18,17) = 1;
graph(22,17) = 1;
graph(21,18) = 1;
graph(21,18) = 1;
graph(20,19) = 1;
graph(20,19) = 1;
graph(23,20) = 1;
graph(23,20) = 1;
graph(22,21) = 1;

% Transformer data
graph(3,24) = 1;
graph(10,5) = 1;
graph(9,11) = 1;
graph(12,9) = 1;
graph(10,11) = 1;
graph(10,12) = 1;

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

topology = "ieee";
graph = sparse(graph);
%unknowns_tensor = sparse(unknowns_tensor); Matlab doesnt support sparse 3D
save('ieee24.mat', 'graph', 'unknowns_tensor', 'topology', 'L');
excited = zeros(L,1);
measured = zeros (L,1);

% Display networks with unknowns
if display
    for i = 1:ntrials
        main_identifiable("IEEE24", graph, unknowns_tensor(:,:,i), excited, measured, 1, 0);
    end
end

end