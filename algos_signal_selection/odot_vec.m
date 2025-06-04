% a,b,c are vectors of same size which elements belong to {-1, 0, 1}

function c = odot_vec(a,b)
    n = length(a);
    c = zeros(size(a));
    for i = 1:n
        c(i) = odot(a(i), b(i));
    end
end