% Given a graph and its excited and measured sets of nodes, this function
% calls identifiable.m and prints its results on the command window.
% It then also calls plot_identifiable to plot the results.
%
% For generic decoupled identifiability, we also compute and plot results
% for the generic local identifiability on the decoupled network.
%
% For explanations about input/output, see specifications of identifiable.m

function main_identifiable(name, graph, unknowns, excited, measured, nsamples, decoupled)

fprintf('%s ', name);
switch decoupled
    case 1 % Generic decoupled identifiability
        fprintf('(decoupled identifiability): ');
        [net, id_modules, nid_modules] = identifiable(graph, unknowns, excited, measured, nsamples, decoupled);
        plot_identifiable(name, graph, excited, measured, nsamples, decoupled, id_modules, nid_modules);

        if net
            fprintf('Identifiable\n');
        else
            fprintf('Not identifiable\n');
            fprintf('Could not identify %d modules:\n', size(nid_modules,1));
            disp(nid_modules);
        end

        % Generic local identifiability of the decoupled network
        fprintf('%s ', name);
        fprintf('(local identifiability on decoupled network): ');
        [L,~] = size(graph);
        dec_graph = [graph zeros(L,L); unknowns graph];
        dec_unknowns = [zeros(L,L) zeros(L,L); unknowns zeros(L,L)];
        dec_excited = [excited; zeros(L,1)];
        dec_measured = [zeros(L,1); measured];
        [net, id_modules, nid_modules] = identifiable(dec_graph, dec_unknowns,...
            dec_excited, dec_measured, nsamples, 0);
        dec_name = strcat(name, ' - Decoupled network');
        plot_identifiable(dec_name, dec_graph, dec_excited, dec_measured,...
            nsamples, 0, id_modules, nid_modules);

        if net
            fprintf('Identifiable\n');
        else
            fprintf('Not identifiable\n');
            fprintf('Could not identify %d modules:\n', size(nid_modules,1));
            disp(nid_modules);
        end
     
    otherwise % Generic local identifiability
        fprintf('(local identifiability): ');
        [net, id_modules, nid_modules] = identifiable(graph, unknowns, excited, measured, nsamples, decoupled);
        plot_identifiable(name, graph, excited, measured, nsamples, decoupled, id_modules, nid_modules);
        if net
            fprintf('Identifiable\n');
        else
            fprintf('Not identifiable\n');
            fprintf('Could not identify %d modules:\n', size(nid_modules,1));
            disp(nid_modules);
        end
end

end