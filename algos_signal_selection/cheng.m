% Algo of pseudotree merging allocation
% Cf paper : Xiaodong Cheng, Shengling Shi, and Paul MJ Van den Hof.
% Allocation of excitation signals for generic identifiability of
% linear dynamic networks. IEEE Transactions on Automatic Control

function excited = cheng(graph)

L = length(graph);

% Permutation such that all zero columns are at the end
zeroCols = all(graph == 0, 1);  % Logical array for zero columns
nonZeroCols = ~zeroCols;
% Rearrange columns and rows to move zero columns to the end
permutation = [find(nonZeroCols), find(zeroCols)];
graph_permuted = graph(permutation, permutation);

% Initialization of covering with minimal pseudotrees
sinks = find(all(graph_permuted == 0, 1));
roots_init = 1:L;
roots_init(ismember(roots_init, sinks)) = [];
n = numel(roots_init);
covering = zeros(L,L,n);
for i = 1:n
    covering(:,roots_init(i),i) = graph_permuted(:,roots_init(i));
end

% Construction of characteristic matrix
M = zeros(n,n);
graph_complex = graph_permuted + 1i * eye(L);
for i = 1:n
    for j = 1:n
        a = graph_complex(:,i).' * graph_complex(:,j);
        if i == j
            M(i,j) = 0;
        elseif real(a) == 0 && imag(a) ~= 0 && graph_permuted(i,j) ~= 0
            M(i,j) = 1;
        elseif a == 0
            M(i,j) = -1;
        else %if real(a) ~= 0 || (real(a) == 0 && imag(a) ~= 0 && graph_permuted(i,j) == 0)
            M(i,j) = 0;
        end
    end
end

%% Algorithm 1 : Disjoint Pseudotree Merging

% First loop
nb_lonely1 = 1;
while nb_lonely1 > 0
    % Find matrix entry which is the only 1 in its row
    idx_lonely1 = [];
    for row = 1:size(M, 1)
        % Check if the row contains exactly one '1'
        if sum(M(row, :) == 1) == 1
            % Find the column index of the '1' in the row
            col = find(M(row, :) == 1);
            % Store the (i,j) pair
            idx_lonely1 = [idx_lonely1; row, col];
        end
    end
    nb_lonely1 = size(idx_lonely1, 1);

    if nb_lonely1 > 0
        i_new = idx_lonely1(1,1);
        j_new = idx_lonely1(1,2);
        % If several ones, take the one with most -1 entries in its row
        if nb_lonely1 > 1
            max_minus1 = -1;
            for idx = nb_lonely1
                row = idx_lonely1(idx,1);
                col = idx_lonely1(idx,2);
                nb_minus1 = sum(M(row,:) == -1);
                if nb_minus1 > max_minus1
                    max_minus1 = nb_minus1;
                    i_new = row;
                    j_new = col;
                end
            end
        end
        % Update M and covering
        [covering, M] = merge_covering(covering, M, i_new, j_new);
    end
end

% Second loop
mergeable = true;
while mergeable
    max_minus1 = -1;
    mergeable = false;
    for row = 1:size(M, 1)
        idx_1 = find(M(row, :) == 1);
        if numel(idx_1) > 0
            mergeable = true;
            nb_minus1 = sum(M(row,:) == -1);
            if nb_minus1 > max_minus1
                max_minus1 = nb_minus1;
                i_new = row;
                j_new = idx_1(1);
            end
        end
    end
    if mergeable
        % Update M and covering
        [covering, M] = merge_covering(covering, M, i_new, j_new);
    end
end

% Print covering for debugging
% ncov = size(covering,3);
% for i = 1:ncov
%     main_identifiable("Pseudotree " + i, covering(:,:,i), zeros(L,L), zeros(L,1), zeros(L,1), 1, 0);
% end

%% Algorithm 2 : Allocation of Excitation Signals

% Select a root for each pseudotree
% Definition: v is a root if there is exactly 1 path from v to every other 
% vertex of the pseudotree
ncov = size(covering,3);
roots = zeros(ncov,1);
for i = 1:ncov
    rowsWithNonZeros = any(covering(:,:,i) ~= 0, 2);
    colsWithNonZeros = any(covering(:,:,i) ~= 0, 1);
    vertices = find(rowsWithNonZeros | colsWithNonZeros');
    nptree = numel(vertices);
    T = inv(eye(nptree) - rand(nptree,nptree) .* covering(vertices,vertices,i));
    for j = 1:nptree
        % If j is a root of pseudotree i
        if all(T(:,j) ~= 0)
            roots(i) = vertices(j);
        end
    end
end

% Deprecated version in comments: error from Cheng paper:
% Not only the nodes of each pseudotree must be considered in the vdp if
% condition, but all the nodes of the graph. Indeed, removing a root
% influences other nodes than just its outneighbours (or 1st generation
% children): it influences all its grandchildren.

% Remove unnecessary roots
T = inv(eye(L) - rand(L,L) .* graph_permuted);
nb_removed = 0;
% Loop on the pseudotrees
for i = 1:ncov
    roots_new = roots;
    roots_new(i - nb_removed) = [];
    % rowsWithNonZeros = any(covering(:,:,i) ~= 0, 2); % depr
    % colsWithNonZeros = any(covering(:,:,i) ~= 0, 1); % depr
    % vertices = find(rowsWithNonZeros | colsWithNonZeros'); % depr
    % nptree = numel(vertices); % depr
    remove = true;
    % Loop on the vertices of pseudotree i
    for j = 1:L % 1:nptree % depr
        % inneighbours = find(graph_permuted(vertices(j),:)); % depr
        inneighbours = find(graph_permuted(j,:));
        b = rank(T(inneighbours, roots_new)); % generic rank
        if b ~= numel(inneighbours)
            remove = false;
            break;
        end
    end
    if remove
        %disp("===============================");
        %disp("Pseudotree number " + i);
        %disp("We remove " + roots(i - nb_removed));
        roots = roots_new;
        nb_removed = nb_removed + 1;
    end
end

excited_permuted = zeros(L,1);
excited_permuted(roots) = 1;

% Make inverse permutation
inverse_permutation = zeros(L,1);
inverse_permutation(permutation) = 1:L;
excited = excited_permuted(inverse_permutation);

%nB_cheng = sum(excited);
%main_identifiable("Cheng's solution permuted", graph_permuted, graph_permuted, excited_permuted, (1:L)', 1, 0);
%main_identifiable("Cheng's solution", graph, graph, excited, (1:L)', 1, 0);

end

