% Given a graph and its excited/measured nodes, calls identifiable and
% prints its results on the command window. Also calls plot_identifiable to
% plot the results.

% For explanations about input/output, see specifications of identifiable.m

function main_identifiable(name, graph, excited, measured, nsamples)

fprintf('%s: ', name);

[net, i_edges, ni_edges] = identifiable(graph, excited, measured, nsamples);
plot_identifiable(name, graph, excited, measured, nsamples, i_edges, ni_edges);

if net
    fprintf('Identifiable\n');
else
    fprintf('Not identifiable\n');
    fprintf('Could identify %d edges:\n', size(i_edges,1));
    disp(i_edges);
    fprintf('Could not identify %d edges:\n', size(ni_edges,1));
    disp(ni_edges);
end

end