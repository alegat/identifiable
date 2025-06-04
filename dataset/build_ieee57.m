% Encodes a Standard IEEE Power Flow Test Cases
% cfr https://icseg.iti.illinois.edu/power-cases/
% Makes 10 selections of 1/3 of the edges as unknown modules.
% Stores the networks in a .m file
% If display == 1, diplays each network generated.
% The random seed is fixed, so each call to this function will output the
% same networks.

function build_ieee57(display)
close all;

ratio = 1/3;
ntrials = 10;
L = 57;
graph = zeros(L, L);

% Branch data
graph(2,1) = 1;
graph(15,1) = 1;
graph(16,1) = 1;
graph(17,1) = 1;
graph(3,2) = 1;
graph(4,3) = 1;
graph(15,3) = 1;
graph(5,4) = 1;
graph(6,4) = 1;
graph(6,5) = 1;
graph(7,6) = 1;
graph(8,6) = 1;
graph(8,7) = 1;
graph(9,8) = 1;
graph(10,9) = 1;
graph(11,9) = 1;
graph(12,9) = 1;
graph(13,9) = 1;
graph(12,10) = 1;
graph(13,11) = 1;
graph(13,12) = 1;
graph(16,12) = 1;
graph(17,12) = 1;
graph(14,13) = 1;
graph(15,13) = 1;
graph(15,14) = 1;
graph(19,18) = 1;
graph(20,19) = 1;
graph(22,21) = 1;
graph(23,22) = 1;
graph(38,22) = 1;
graph(24,23) = 1;
graph(25,24) = 1;
graph(25,24) = 1;
graph(30,25) = 1;
graph(27,26) = 1;
graph(28,27) = 1;
graph(29,28) = 1;
graph(52,29) = 1;
graph(31,30) = 1;
graph(32,31) = 1;
graph(33,32) = 1;
graph(35,34) = 1;
graph(36,35) = 1;
graph(37,36) = 1;
graph(40,36) = 1;
graph(38,37) = 1;
graph(39,37) = 1;
graph(44,38) = 1;
graph(48,38) = 1;
graph(49,38) = 1;
graph(42,41) = 1;
graph(43,41) = 1;
graph(41,56) = 1;
graph(42,56) = 1;
graph(45,44) = 1;
graph(47,46) = 1;
graph(48,47) = 1;
graph(49,48) = 1;
graph(50,49) = 1;
graph(51,50) = 1;
graph(53,52) = 1;
graph(54,53) = 1;
graph(55,54) = 1;
graph(56,57) = 1;

% Transformer data
graph(18,4) = 1;
graph(18,4) = 1;
graph(29,7) = 1;
graph(55,9) = 1;
graph(51,10) = 1;
graph(41,11) = 1;
graph(43,11) = 1;
graph(49,13) = 1;
graph(46,14) = 1;
graph(45,15) = 1;
graph(20,21) = 1;
graph(26,24) = 1;
graph(32,34) = 1;
graph(57,39) = 1;
graph(56,40) = 1;

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
save('ieee57.mat', 'graph', 'unknowns_tensor', 'topology', 'L');
excited = zeros(L,1);
measured = zeros (L,1);

% Display networks with unknowns
if display
    for i = 1:ntrials
        main_identifiable("IEEE57", graph, unknowns_tensor(:,:,i), excited, measured, 1, 0);
    end
end

end