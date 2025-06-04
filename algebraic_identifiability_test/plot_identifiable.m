% Plots the graph, highlights the excited and measured nodes, as well as the
% identifiable and non-identifiable modules given as argument.
%
% Legend:
% Nodes are marked as:
%   v (yellow) if excited
%   ^ (orange) if measured
%   * (brown) if excited and measured (superposition of v and ^)
%   o if neither excited nor measured
% Edges are coloured as:
%   green if identifiable
%   red if non-identifiable
%   blue if known, or not specified (nargin <= 6)

function plot_identifiable(name, graph, excited, measured, nsamples, ...
    relax, id_modules, nid_modules)

[t, s] = find(graph);
G = digraph(s, t);
fontsize = 16;

% Colours
green = [0.4660 0.6740 0.1880];
yellow = [0.9290 0.6940 0.1250];
orange = [0.8500 0.3250 0.0980];
brown = [0.6350 0.0780 0.1840];
blue = [0 0.4470 0.7410];

figure;
p = plot(G);
p.Marker = 'o';
p.NodeColor = blue;
p.MarkerSize = 13;
p.LineWidth = 3;
p.ArrowSize = 13;
p.NodeFontSize = fontsize;

% Highlight nodes
highlight(p, find(measured), 'Marker', '^', 'NodeColor', orange);
highlight(p, find(excited), 'Marker', 'v', 'NodeColor', yellow);
highlight(p, find(excited .* measured),'Marker','h','MarkerSize',p.MarkerSize+5 ...
          , 'NodeColor', brown);

if nargin > 6
    % Highlight identifiable modules
    if ~isscalar(id_modules)
        for i = 1:size(id_modules,1)
            highlight(p, [id_modules(i,2) id_modules(i,1)], 'EdgeColor', green);
        end
    end
    % Highlight non-identifiable modules
    if ~isscalar(nid_modules)
        for i = 1:size(nid_modules,1)
            highlight(p, [nid_modules(i,2) nid_modules(i,1)], 'EdgeColor', 'red');
        end
    end
end

% Title
switch relax
    case 1
        relaxName = " - Generic decoupled identifiability ";
    otherwise
        relaxName = " - Generic local identifiability ";
end
if nsamples == 0
    nsamplesName = "- Symbolic";
elseif nsamples == -1
    nsamplesName = "";
else
    nsamplesName = sprintf("- %d samples", nsamples);
end
title(name + relaxName + nsamplesName,'FontSize',fontsize);


% Legend
hold on;
leg = plot(NaN,NaN,'v',NaN,NaN,'^',NaN,NaN,'h',NaN,NaN,'o',NaN,NaN,'g',NaN,NaN,'r',NaN,NaN);
leg(1).MarkerSize = 10;
leg(1).Color = yellow;
leg(1).MarkerFaceColor = yellow;
leg(2).MarkerSize = 10;
leg(2).Color = orange;
leg(2).MarkerFaceColor = orange;
leg(3).MarkerSize = 15;
leg(3).Color = brown;
leg(3).MarkerFaceColor = brown;
leg(4).MarkerSize = 10;
leg(4).Color = blue;
leg(4).MarkerFaceColor = blue;
leg(5).Color = green;
leg(5).LineWidth = 2;
leg(6).LineWidth = 2;
leg(7).LineWidth = 2;
leg(7).Color = blue;
legend(leg,'Excited','Measured','Excited & measured','Passive',...
    'Identifiable','Non-identifiable','Known','FontSize',fontsize,...
    'Location','southeast');
%legend('boxoff')

end 