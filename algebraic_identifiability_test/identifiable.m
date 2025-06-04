% Returns the generic local identifiability of each unknown module of a graph,
% given its excited and measured nodes.

% INPUT
% graph: directed adjacency matrix. Size LxL
% entry (i,j) = 1 if there is an edge from j to i
%             = 0 otherwise
% unknowns: directed adjacency matrix of unknown modules. Size LxL
% entry (i,j) = 1 if there is an *unknown* edge from j to i
%             = 0 otherwise
% excited: binary selection of excited nodes. Size Lx1
%     entry = 1 if excited
%           = 0 otherwise
% measured: binary selection of measured nodes. Size Lx1
%     entry = 1 if measured
%           = 0 otherwise
% nsamples: number of random samples.
%           = 0 for symbolic computation
% decoupled = 1 for generic decoupled identifiability
%         otherwise, it is generic local identifiability
%
% OUTPUT
% net = 1 if all unknown modules are generically locally identifiable
%     = 0 otherwise
% i_modules is a m x 2 array where m is the nb of generically locally identifiable modules
%            (i,1) is the outgoing node of edge i
%            (i,2) is the ingoing node of edge i
% ni_modules = 0 if all modules are generically locally identifiable
%            otherwise, it is a n x 2 array where n is the nb of non-identifiable modules
%            (i,1) is the outgoing node of edge i
%            (i,2) is the ingoing node of edge i

function [net, id_modules, nid_modules] = identifiable(graph, unknowns, ...
    excited, measured, nsamples, decoupled)

switch nargin
    case 3
        nsamples = 0;
        decoupled = 0;
    case 4
        decoupled = 0;
end

tol = 1e-9;
[L,~] = size(graph);
B = diag(excited); [~,Bcol] = find(B); B = B(:,Bcol);
C = diag(measured);[Crow,~] = find(C); C = C(Crow,:);
Id = eye(L);
[unknowns_row, unknowns_col] = find(unknowns);
nbunknowns = sum(unknowns(:));

% Initialization: 1 for identifiable; 0 for non-identifiable
net = 0;
modules = zeros(nbunknowns,1);

if nsamples > 0
    nsamples = nsamples - 1;
    symbolic = 0;
else
    symbolic = 1;
end

for j = 0:nsamples
    if symbolic
        syms G [L L]
    else
        G = rand(L,L);
    end
    G = G .* graph;
    if decoupled
        if symbolic
            syms Gbis [L L]
        else
            Gbis = rand(L,L);
        end
        Gbis = Gbis .* graph;
        K = kron_submatrix(((Id-G)\B)', C/(Id-Gbis), unknowns);
    else
        K = kron_submatrix(((Id-G)\B)', C/(Id-G), unknowns);
    end

    r = rank(K);
    if r == nbunknowns
        net = 1;  % if the network is locally identifiable for 1 sample,
        break;    % then it is generically locally identifiable
    else
        V = null(K); % orthonormal basis for the null space of K
                     % 'zeros' mostly < 1e-15, 'nonzeros' mostly > 1e-3
        if symbolic
            Vlogic = V ~= 0;
        else
            Vlogic = abs(V) > tol; % cast numerical array to logical binary
        end
        v = ~any(Vlogic,2); % 1 if identifiable; 0 if non-identifiable
        modules = modules + v/(nsamples+1);
    end
end
modules = round(modules) ~= 0;

% Computation of identifiable and non-identifiable unknown modules
if net
    id_modules = [unknowns_row unknowns_col];
    nid_modules = 0;
else
    id_row = unknowns_row(modules);
    id_col = unknowns_col(modules);
    
    nid_row = unknowns_row(~modules);
    nid_col = unknowns_col(~modules);
    
    id_modules = [id_row id_col]; % identifiable modules
    nid_modules = [nid_row nid_col]; % non-identifiable modules
end

end