% Encodes a Standard IEEE Power Flow Test Cases
% cfr https://icseg.iti.illinois.edu/power-cases/
% Makes 10 selections of 1/3 of the edges as unknown modules.
% Stores the networks in a .m file
% If display == 1, diplays each network generated.
% The random seed is fixed, so each call to this function will output the
% same networks.

function build_ieee30(display)
close all;

ratio = 1/3;
ntrials = 10;
L = 30;
graph = zeros(L, L);

% Branch data
graph(2,1) = 1;
graph(3,1) = 1;
graph(4,2) = 1;
graph(5,2) = 1;
graph(6,2) = 1;
graph(4,3) = 1;
graph(6,4) = 1;
graph(7,5) = 1;
graph(7,6) = 1;
graph(8,6) = 1;
graph(28,6) = 1;
graph(28,8) = 1;
graph(10,9) = 1;
graph(11,9) = 1;
graph(17,10) = 1;
graph(20,10) = 1;
graph(21,10) = 1;
graph(22,10) = 1;
graph(13,12) = 1;
graph(14,12) = 1;
graph(15,12) = 1;
graph(16,12) = 1;
graph(15,14) = 1;
graph(18,15) = 1;
graph(23,15) = 1;
graph(17,16) = 1;
graph(19,18) = 1;
graph(20,19) = 1;
graph(22,21) = 1;
graph(24,22) = 1;
graph(24,23) = 1;
graph(25,24) = 1;
graph(26,25) = 1;
graph(27,25) = 1;
graph(29,27) = 1;
graph(30,27) = 1;
graph(30,29) = 1;

% Transformer data
graph(12,4) = 1;
graph(9,6) = 1;
graph(10,6) = 1;
graph(27,28) = 1;

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
save('ieee30.mat', 'graph', 'unknowns_tensor', 'topology', 'L');
excited = zeros(L,1);
measured = zeros (L,1);

% Display networks with unknowns
if display
    for i = 1:ntrials
        main_identifiable("IEEE30", graph, unknowns_tensor(:,:,i), excited, measured, 1, 0);
    end
end

end