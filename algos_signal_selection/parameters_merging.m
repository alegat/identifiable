function [nrun, stop] = parameters_merging(L)

nrun = 10;

% Watts style
if L <= 50
    stop = L/50;
else
    stop = 1;
end

end