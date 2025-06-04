% Computes kron(BTt, CT) * I° efficiently
% i.e. only the columns of kron(BTt, CT) corresponding to unknowns.
%
% Faster and less memory usage than computing the full kronecker product,
% and then only selecting the columns corresponding to unknowns.
%
% Time & memory : O(nB * nC * m°) instead of O(nB * nC * n^2).

function K = kron_submatrix(BTt, CT, unknowns)

[nB, ~] = size(BTt);
[nC, ~] = size(CT);
nbunknowns = sum(unknowns, 'all');
[unknowns_row, unknowns_col] = find(unknowns);
K = zeros(nB * nC, nbunknowns);
if isa(CT, 'sym')
    K = sym(K);
end

for b = 1:nB
    for c = 1:nC
        for l = 1:nbunknowns
            start_l = unknowns_col(l);
            end_l = unknowns_row(l);
            K((b-1)*nC + c,l) = BTt(b, start_l) * CT(c, end_l);
        end
    end
end

end