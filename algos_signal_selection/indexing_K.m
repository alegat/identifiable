% Computes the row indices corresponding to the submatrix kron(B'T',CT) of
% kron(T',T). In other words, those are the indices corresponding to 
% excitations excited_idx and measurements measured_idx in Kfull.

function indices = indexing_K(L, excited_idx, measured_idx)

nB = numel(excited_idx);
nC = numel(measured_idx);

indices = zeros(nB*nC,1);
for i = 1:nB*nC
    indices(i) = (excited_idx(ceil(i/nC)) - 1) * L ...
                + measured_idx(mod(i-1,nC) + 1);
end

end