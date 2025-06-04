function [covering, M_new] = merge_covering(covering, M, i_new, j_new)

% Update M
M_new = M;
M_new(j_new,:) = odot_vec(M(i_new,:), M(j_new,:));
M_new(:,j_new) = odot_vec(M(:,i_new), M(:,j_new));
M_new(i_new,:) = [];
M_new(:,i_new) = [];

% Update covering
covering(:,:,j_new) = covering(:,:,i_new) + covering(:,:,j_new);
covering(:,:,i_new) = [];

end