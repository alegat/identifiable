% Computes the identifiability of a graph given its excited and measured
% nodes. Returns the local generic identifiability of each edge.

% INPUT
% graph: directed adjacency matrix. Size LxL
%           = 1 if there is an edge from j to i
%           = 0 otherwise
% excited: binary selection of excited nodes. Size Lx1
%           = 1 if excited
%           = 0 otherwise
% measured: binary selection of measured nodes. Size Lx1
%           = 1 if measured
%           = 0 otherwise
% nsamples: number of random samples.
%           = 0 for symbolic computation
% OUTPUT
% net = 1 if all edges are generically identifiable
%     = 0 otherwise
% i_edges is a m x 2 array where m is the nb of identifiable edges.
%            (i,1) is the outgoing node of identifiable edge i
%            (i,2) is the ingoing node of identifiable edge i
% ni_edges = 0 if all edges are generically identifiable
%            otherwise is a m x 2 array where m is the nb of
%            non-identifiable edges.
%            (i,1) is the outgoing node of n-i edge i
%            (i,2) is the ingoing node of n-i edge i

function [net, i_edges, ni_edges] = identifiable(graph, excited, measured,...
    nsamples)

switch nargin
    case 3
        nsamples = 0;
end

tol = 1e-9;
[L,~] = size(graph); % nb nodes
B = diag(excited); [~,Bcol] = find(B); B = B(:,Bcol);
C = diag(measured);[Crow,~] = find(C); C = C(Crow,:);
Id = eye(L);
delta = graph(:)';
deltaLogic = delta ~= 0; % logical version of delta
[edges_row, edges_col] = find(graph);
nbedges = sum(delta);

% Initialization: 1 if identifiable; 0 if non-identifiable
net = 0;
edges = zeros(nbedges,1);

if nsamples > 0
    nsamples = nsamples - 1;
    symbolic = 0;
else
    symbolic = 1;
end
for i = 0:nsamples
    if symbolic
        syms A [L L]
    else
        A = rand(L,L);
    end
    G = A .* graph;
    K0 = kron(((Id-G)\B)',C/(Id-G));
    K = K0(:,deltaLogic); % only keep the columns matching nonzero delta(i)
    r = rank(K);
    if r == nbedges
        net = 1; % if network identif for 1 sample, then generic identif
    else
        N = null(K); % orthonormal basis for the null space of Kr
                      % 'zeros' mostly <1e-15, 'nonzeros' mostly >1e-3
        if symbolic
            Nlogic = N ~= 0;
        else
            Nlogic = abs(N) > tol; % cast numerical array to logical binary
        end
        n = ~any(Nlogic,2); % 1 if identifiable; 0 if non-identifiable
        edges = edges | n; % cumulated n. If one edge is identifiable for
                           % one sample, it is generically identifiable
    end
end

% Computation of identifiable and non-identifiable edges
if net
    i_edges = [edges_row edges_col];
    ni_edges = 0;
else
    i_row = edges_row(edges);
    i_col = edges_col(edges);
    
    ni_row = edges_row(~edges);
    ni_col = edges_col(~edges);
    
    i_edges = [i_row i_col]; % identifiable edges
    ni_edges = [ni_row ni_col]; % non-identifiable edges
end

end