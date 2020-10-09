% Plots the graph, excited/measured nodes and (non)-identifiable edges
% given as argument.

% Nodes are marked as:
%   v if excited
%   ^ if measured
%   a star if excited and measured (superposition of v and ^)
%   o if neither excited nor measured
% Edges are coloured as:
%   green if identifiable
%   red if non-identifiable
%   blue if not specified (nargin <= 6)

function plot_identifiable(name, graph, excited, measured, nsamples, ...
      i_edges, ni_edges)

[t, s] = find(graph);
G = digraph(s, t);
fontsize = 16;
myBlue = [0 0.4470 0.7410];
myGreen = [0.4660 0.6740 0.1880];

figure;
p = plot(G);
p.Marker = 'o';
p.NodeColor = myBlue;
p.MarkerSize = 10;
p.LineWidth = 2;
p.ArrowSize = 10;
p.NodeFontSize = fontsize; % comment this line if old version of Matlab

% Highlight nodes
highlight(p, find(measured), 'Marker', '^');
highlight(p, find(excited), 'Marker', 'v');
highlight(p, find(excited .* measured),'Marker','h','MarkerSize',15);

if nargin > 6
    % Highlight identifiable edges
    if ~isscalar(i_edges)
        for i = 1:size(i_edges,1)
            highlight(p, [i_edges(i,2) i_edges(i,1)], 'EdgeColor',myGreen);
        end
    end
    % Highlight non-identifiable edges
    if ~isscalar(ni_edges)
        for i = 1:size(ni_edges,1)
            highlight(p, [ni_edges(i,2) ni_edges(i,1)], 'EdgeColor', 'red');
        end
    end
end

% Title
if nsamples == 0
    nsamplesName = '- Symbolic';
else
    nsamplesName = sprintf('- %d samples', nsamples);
end
title([name, ' ', nsamplesName],'FontSize',fontsize);

% Legend
hold on;
leg = plot(NaN,NaN,'v',NaN,NaN,'^',NaN,NaN,'h',NaN,NaN,'g',NaN,NaN,'r');
for i = 1:3
    leg(i).MarkerSize = 10;
    leg(i).Color = myBlue;
    leg(i).MarkerFaceColor = myBlue;
end
leg(3).MarkerSize = 15;
leg(4).Color = myGreen;
leg(4).LineWidth = 2;
leg(5).LineWidth = 2;
legend(leg,{'Excited','Measured','Excited & measured',...
    'Identifiable','Non-identifiable'},'FontSize',fontsize);

end 