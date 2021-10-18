% Trying to find networks which are decoupled-identifiable but not locally
% identifiable.
%
% In main_identifiable(name, graph, excited, measured, nsamples, relax):
% - nsamples = 10 is good, 20 is better to avoid numerical errors due to
% isolated nodes and loops
% - relax = 0 for local identifiability
% - relax = 1 for decoupled identifiability
%
% Careful : this code tends to see fake counterexamples when there are isolated
% nodes or loops in the network, because of bad numerical behaviour.
%
% If the code does not display anything, no counterexample has been found

close all;

% n = nb of nodes
n = 10; % between 10 and 100 (30 is good for 1e3 networks)

% coef_G = density of edges
coef_G = 0.6; % between 0.6 and 1 (0.6 is good)

% coef_E = density of excitations and measures
coef_E = 2; % between 1.5 and 4 (2 is good)

% nbnet = number of randomly generated networks
nbnet = 1e4;

% Recommended values for a reasonable running time
% n = 3 for 1e6 networks
% n = 10 for 1e5 networks
% n = 20 for 1e4 networks
% n = 30 for 1e3 networks
% n = 100 for 2 networks

tic
for i = 1:nbnet
    graph = round(coef_G*rand(n,n));
    measured = min(1,round(coef_E*rand(n,1)));
    excited = min(1,round(coef_E*rand(n,1)));
    [~, i_edges0, ~] = identifiable(graph, excited, measured, 20, 0);
    [~, i_edges1, ~] = identifiable(graph, excited, measured, 20, 1);
    if length(i_edges0) ~= length(i_edges1)
        disp("==================================================");
        disp("============= COUNTER-EXAMPLE FOUND ==============");
        disp("Network decoupled-identif but not locally identif:");
        disp("==================================================");
        disp(measured);
        disp(excited);
        disp(i_edges0);
        disp(i_edges1);
        % 100 samples are taken to make sure to avoid numerical errors
        main_identifiable('Test', graph, excited, measured, 100, 0);
        main_identifiable('Test', graph, excited, measured, 100, 1);
    end
end
toc
