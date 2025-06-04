% Encodes a Standard IEEE Power Flow Test Cases
% cfr https://icseg.iti.illinois.edu/power-cases/
% Makes 10 selections of 1/3 of the edges as unknown modules.
% Stores the networks in a .m file
% If display == 1, diplays each network generated.
% The random seed is fixed, so each call to this function will output the
% same networks.

function build_ieee118(display)
close all;

ratio = 1/3;
ntrials = 10;
L = 118;
graph = zeros(L,L);

% Branch data
graph(2,1) = 1;
graph(3,1) = 1;
graph(12,2) = 1;
graph(5,3) = 1;
graph(12,3) = 1;
graph(5,4) = 1;
graph(11,4) = 1;
graph(6,5) = 1;
graph(11,5) = 1;
graph(7,6) = 1;
graph(12,7) = 1;
graph(9,8) = 1;
graph(30,8) = 1;
graph(10,9) = 1;
graph(12,11) = 1;
graph(13,11) = 1;
graph(14,12) = 1;
graph(16,12) = 1;
graph(117,12) = 1;
graph(15,13) = 1;
graph(15,14) = 1;
graph(17,15) = 1;
graph(19,15) = 1;
graph(33,15) = 1;
graph(17,16) = 1;
graph(18,17) = 1;
graph(31,17) = 1;
graph(113,17) = 1;
graph(19,18) = 1;
graph(20,19) = 1;
graph(34,19) = 1;
graph(21,20) = 1;
graph(22,21) = 1;
graph(23,22) = 1;
graph(24,23) = 1;
graph(25,23) = 1;
graph(32,23) = 1;
graph(70,24) = 1;
graph(72,24) = 1;
graph(27,25) = 1;
graph(30,26) = 1;
graph(28,27) = 1;
graph(32,27) = 1;
graph(115,27) = 1;
graph(29,28) = 1;
graph(31,29) = 1;
graph(38,30) = 1;
graph(32,31) = 1;
graph(113,32) = 1;
graph(114,32) = 1;
graph(37,33) = 1;
graph(36,34) = 1;
graph(37,34) = 1;
graph(43,34) = 1;
graph(36,35) = 1;
graph(37,35) = 1;
graph(39,37) = 1;
graph(40,37) = 1;
graph(65,38) = 1;
graph(40,39) = 1;
graph(41,40) = 1;
graph(42,40) = 1;
graph(42,41) = 1;
graph(49,42) = 1;
graph(49,42) = 1;
graph(44,43) = 1;
graph(45,44) = 1;
graph(46,45) = 1;
graph(49,45) = 1;
graph(47,46) = 1;
graph(48,46) = 1;
graph(49,47) = 1;
graph(69,47) = 1;
graph(49,48) = 1;
graph(50,49) = 1;
graph(51,49) = 1;
graph(54,49) = 1;
graph(54,49) = 1;
graph(66,49) = 1;
graph(66,49) = 1;
graph(69,49) = 1;
graph(57,50) = 1;
graph(52,51) = 1;
graph(58,51) = 1;
graph(53,52) = 1;
graph(54,53) = 1;
graph(55,54) = 1;
graph(56,54) = 1;
graph(59,54) = 1;
graph(56,55) = 1;
graph(59,55) = 1;
graph(57,56) = 1;
graph(58,56) = 1;
graph(59,56) = 1;
graph(59,56) = 1;
graph(60,59) = 1;
graph(61,59) = 1;
graph(61,60) = 1;
graph(62,60) = 1;
graph(62,61) = 1;
graph(66,62) = 1;
graph(67,62) = 1;
graph(64,63) = 1;
graph(65,64) = 1;
graph(68,65) = 1;
graph(67,66) = 1;
graph(81,68) = 1;
graph(116,68) = 1;
graph(70,69) = 1;
graph(75,69) = 1;
graph(77,69) = 1;
graph(71,70) = 1;
graph(74,70) = 1;
graph(75,70) = 1;
graph(72,71) = 1;
graph(73,71) = 1;
graph(75,74) = 1;
graph(77,75) = 1;
graph(118,75) = 1;
graph(77,76) = 1;
graph(118,76) = 1;
graph(78,77) = 1;
graph(80,77) = 1;
graph(80,77) = 1;
graph(82,77) = 1;
graph(79,78) = 1;
graph(80,79) = 1;
graph(96,80) = 1;
graph(97,80) = 1;
graph(98,80) = 1;
graph(99,80) = 1;
graph(83,82) = 1;
graph(96,82) = 1;
graph(84,83) = 1;
graph(85,83) = 1;
graph(85,84) = 1;
graph(86,85) = 1;
graph(88,85) = 1;
graph(89,85) = 1;
graph(87,86) = 1;
graph(89,88) = 1;
graph(90,89) = 1;
graph(90,89) = 1;
graph(92,89) = 1;
graph(92,89) = 1;
graph(91,90) = 1;
graph(92,91) = 1;
graph(93,92) = 1;
graph(94,92) = 1;
graph(100,92) = 1;
graph(102,92) = 1;
graph(94,93) = 1;
graph(95,94) = 1;
graph(96,94) = 1;
graph(100,94) = 1;
graph(96,95) = 1;
graph(97,96) = 1;
graph(100,98) = 1;
graph(100,99) = 1;
graph(101,100) = 1;
graph(103,100) = 1;
graph(104,100) = 1;
graph(106,100) = 1;
graph(102,101) = 1;
graph(104,103) = 1;
graph(105,103) = 1;
graph(110,103) = 1;
graph(105,104) = 1;
graph(106,105) = 1;
graph(107,105) = 1;
graph(108,105) = 1;
graph(107,106) = 1;
graph(109,108) = 1;
graph(110,109) = 1;
graph(111,110) = 1;
graph(112,110) = 1;
graph(115,114) = 1;

% Transformer data
graph(5,8) = 1;
graph(17,30) = 1;
graph(25,26) = 1;
graph(37,38) = 1;
graph(59,63) = 1;
graph(61,64) = 1;
graph(66,65) = 1;
graph(69,68) = 1;
graph(80,81) = 1;

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
save('ieee118.mat', 'graph', 'unknowns_tensor', 'topology', 'L');
excited = zeros(L,1);
measured = zeros (L,1);

% Display networks with unknowns
if display
    for i = 1:ntrials
        main_identifiable("IEEE118", graph, unknowns_tensor(:,:,i), excited, measured, 1, 0);
    end
end

end