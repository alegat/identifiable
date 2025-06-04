% a,b,c are vectors of same size which elements belong to {-1, 0, 1}

function c = otimes_vec(a,b)
    n = length(a);
    c = zeros(size(a));
    for i = 1:n
        c(i) = otimes(a(i), b(i));
    end
end