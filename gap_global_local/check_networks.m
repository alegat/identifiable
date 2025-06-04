% Exhaustive search of all networks with L nodes, nb_modules modules and all
% patterns of excitations and measurements such that :
%           nb_modules = nb_excitations * nb_measurements,
% except the fully excited or measured ones.
% We also ignore the disconnected graphs since they can be checked by
% their respective connected components.
%
% For each network, we compute generic local identifiability. If it is 
% generically locally identifiable, we check genericglobal identifiability 
% using the function solve.
%
% S.conditions is a vector of conditions for each solution of S. The
% condition themselves do not matter (sometimes, they are just 'symtrue'),
% but S.conditions gives the number of different discrete solution.
%
% For example, when solving x^2 == 4, S.x is length 2 : [-2, +2].
% S.conditions is also length 2 : [symtrue, symtrue].
% In this example, we could take either length(S.x) or length(S.conditions)
% to detect that we have 2 discrete solutions.
%
% In our case though, we do not know the name of the variables since it
% depends on the topology. But we know that there will always be a condition
% for each discrete solution.
% So length(S.conditions) > 1 is our detection criterion for several
% solutions "discretely" separated.
%
% Typical solution :
% S =  
%           G2_1: z1
%           G3_1: z
%           G3_2: z2
%           G2_3: z3
%        Gbis2_1: z1
%        Gbis3_1: z
%        Gbis3_2: z2
%        Gbis2_3: z3
%     parameters: [z    z1    z2    z3]
%     conditions: z2^2*z3^2 + 1 ~= 2*z2*z3 & z ~= 0 & z1 ~= 0 & z2 …
%
% S.conditions actually gives the lower-dimensional set on which global
% identifiability does not hold.
%
% We have a potential counter-example if :
% - length(S.conditions) > 1 : we have several solutions "discretely"
% separated
% - length(S.parameters) > nb_modules : the solution space has a larger
% dimension than the number of unknowns. Normally it should then not be
% locally identifiable so this case should not arise.
%
% If there is a counter-example, we expect that it will be with
% length(S.parameters) = nb_modules, and it would take the following format
% (e.g. nb_modules = 2) :
% Solution 1 : G = (z, z1), Gbis = (z, z1)
% Solution 2 : G = (z, z1), Gbis = (f(z,z1), g(z,z1))
%                      e.g. Gbis = (2*z + z1, z^2 * z1)
%
% When there is a self-loop in the network we sometimes find a second
% solution which is a particular case of the first one. We don't know why
% this happens but this raises "false" counterexamples. Example of
% solutions we get :
% Solution 1 : G = Gbis = (z, z1)
% Solution 2 : G = Gbis = (z, 1)
% Or 
% Solution 1 : G = Gbis = (z, z1)
% Solution 2 : G = Gbis = (z, 2*z)
% Hence we also exclude self-loops from our search
%
% Graph prompts :
% If L = nb_modules + 1 : Only trees. I would be surprised that a counterex
% arises in a tree topology
% If L > nb_modules + 1 : Disconnected graph. Not interesting to check since
% it can be checked by each of its subnetworks
%
% OUTPUTS :
% check_networks(4,4) : 15 660 networks checked, 1 min elapsed
% check_networks(5,4) : 352 800 networks checked, 6 min elapsed
% nb_modules = 5: either nB=5 and nC=1, or the converse, since nB*nC=nb_modules
% check_networks(< 5,5) : Not possible
% check_networks(5,5) : Would be either full excitation or full measurement
% check_networks(6,5) : 2 488 320 networks checked, 1h37 elapsed
% Still no counter-example found
% check_networks(> 6,5) : Not connected
%
% ===> All networks with <= 5 unknown modules have been checked (except with
% self-loops or nb_modules != nB * nC), and no counter-example has been found
%
% I advise to not go with more modules, it will be untractable

function check_networks(L,nb_modules)

addpath('../algebraic_identifiability_test')

% All graphs of L nodes and nb_modules modules
all_graphs_vec = dec2bin(sum(nchoosek(2.^(0:L^2-1),nb_modules),2)) - '0';
[nb_graphs, ~] = size(all_graphs_vec);

% all_EMP : we will use this matrix for both excitations and measurements
all_EMP = ff2n(L); % EMP : Excitation-measurement pattern
[nb_tot_EMP, ~] = size(all_EMP);
% We remove:
% - first EMP : nothing is excited (resp. measured)
% - last  EMP : full excitation (resp. measurement)
EMPs = all_EMP(2:nb_tot_EMP-1, :);
[nb_EMP, ~] = size(EMPs);

iter = 0;

tic
for i=1:nb_graphs
    graph = reshape(all_graphs_vec(i,:), [L,L]);
    % Loop detection.
    loop = 0;
    for ind = 1:L
        if graph(ind,ind) == 1
            loop = 1;
        end
    end
    % If the graph is not connected, we ignore it.
    % If there is at least one loop, we ignore it.
    if is_connected(graph) && ~loop
        for j=1:nb_EMP
            excited = all_EMP(j,:);
            nb_excitations = sum(excited);
            for k=1:nb_EMP
                measured = all_EMP(k,:);
                nb_measurements = sum(measured);
                if nb_modules == nb_excitations * nb_measurements
                    if mod(iter,1000) == 0
                        disp(iter);
                    end
                    iter = iter+1;
                    % Check generic local identifiability
                    local_bool = identifiable(graph, graph, excited, measured, 5, 0);
                    if local_bool % generically locally identifiable
                        S = solve_CTB_CTbisB(graph, excited, measured);
                        % If length(S.conditions) > 1, it means that there are
                        % multiple discrete solutions
                        potential_counterex = 0;
                        if length(S.conditions) > 1
                            disp("=========================================");
                            disp("======= Potential counter-example =======");
                            disp("======== length(S.conditions) > 1 =======");
                            disp("=========================================");
                            potential_counterex = 1;
                        end
                        % If length(S.parameters) > nb_modules, it means that
                        % it is not generically globally identifiable (but 
                        % it shouldn't be generically locally identifiable neither).
                        if length(S.parameters) > nb_modules
                            disp("=========================================");
                            disp("======= Potential counter-example =======");
                            disp("=== length(S.parameters) > nb_modules ===");
                            disp("=========================================");
                            potential_counterex = 1;
                        end
                        if potential_counterex
                            disp("graph =")
                            disp(graph);
                            fprintf("excited  = ")
                            disp(excited);
                            fprintf("measured = ")
                            disp(measured);
                            % Double check of generic local identifiability
                            % with symbolic computation
                            main_identifiable('test', graph, graph, excited, measured, 0, 0);
                            disp(S);
                        end
                    end
                end
            end
        end
    end
end
fprintf("Nb nodes = %i\n", L);
fprintf("Nb modules = %i\n", nb_modules);
fprintf("Nb of networks checked : %i \n", iter);
toc
end