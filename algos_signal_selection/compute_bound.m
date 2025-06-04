function nbsignals_bound = compute_bound(nbunknowns)

nbsignals_bound = ceil(sqrt(nbunknowns)) ...
                 + ceil(nbunknowns/ceil(sqrt(nbunknowns)));

end