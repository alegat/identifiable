% Builds all networks required by function compare_algos_performance

% Each function call :
% Builds a graph of specific topology
% Makes 10 selections of 1/3 of the edges as unknown modules.
% Stores the networks in a .m file
% If display == 1, diplays each network generated.
% The random seed is fixed, so each call to this function will output the
% same networks.


function build_all
close all;
display = 0;

addpath('../algebraic_identifiability_test')

build_ieee14(display);
build_ieee24(display);
build_ieee30(display);
build_ieee39(display);
build_ieee57(display);
build_ieee118(display);

build_lattice(3, 3, display);
build_lattice(4, 4, display);
build_lattice(4, 5, display);
build_lattice(5, 6, display);
build_lattice(6, 7, display);
build_lattice(7, 7, display);
build_lattice(10, 10, display);

Llist = [10 15 20 30 40 50 100];

for L = Llist
    build_erdos(L, display);
    build_rgg(L, display);
    build_watts(L, display);
end

end