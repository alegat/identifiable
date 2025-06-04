% Merging Measurements - Uniform pass

function [excited_idx, measured_idx_out] = merging_uniform_C(graph, unknowns, ...
    targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
    excited_idx, measured_idx, sigma_B, sigma_C, stop, B_done)

L = length(graph);
nbunknowns = sum(unknowns, 'all');
excited_0 = any(unknowns,1);
measured_0 = any(unknowns,2);
nbsignals_0 = sum(excited_0) + sum(measured_0);

% Stop criterion
% If stop = 0, we do the full merging algo until it is not possible anymore
% If stop = 1, we don't merge at all and go straight to pruning
nbsignals_bound = compute_bound(nbunknowns);
nbsignals_stop = nbsignals_bound + (nbsignals_0 - nbsignals_bound) * stop;
if length(excited_idx) + length(measured_idx) <= nbsignals_stop
    measured_idx_out = measured_idx;
    return
end

% Loop on all measurements, in random order
measured_idx = measured_idx(randperm(length(measured_idx)));
for i = 1:length(measured_idx)
    ci = measured_idx(i);
    outneighbours = find(ancestors(:,ci));
    outneighbours = outneighbours(randperm(length(outneighbours)));

    % Loop on all multi-generational outneighbours of measurement ci
    for j = 1:length(outneighbours)
        nj = outneighbours(j);
        nj_is_measured = ismember(nj, measured_idx);
        inneighbours = find(ancestors(nj,:));
        inneighbours = inneighbours(randperm(length(inneighbours)));
        intersection = intersect(inneighbours, measured_idx);

        % To avoid trying to merge too many measurements that would go
        % below the bound of nbsignals
        nb_merging_max = length(excited_idx) + length(measured_idx) ...
            - nbsignals_bound + 1 - nj_is_measured;
        nb_merging_start = min(nb_merging_max, length(intersection));

        % nb_merging is the nb of excitations we try to merge into nj
        % We are greedy: we start with the highest possible nb_merging and
        % then try with smaller ones
        for nb_merging = nb_merging_start:-1:(2 - nj_is_measured)
            
            if nb_merging == 1
                groups = intersection;
            else
                groups = nchoosek(intersection,nb_merging);
            end
            nb_groups = size(groups, 1);
            groups = groups(randperm(nb_groups),:);

            % Loop on all groups of measurement computed by nchoosek
            for k = 1:nb_groups
                candidates = groups(k,:);
                if nj_is_measured
                    candidates = [nj, candidates];
                end

                % Computation of measured_idx_new
                measured_idx_new = measured_idx;
                if ~ nj_is_measured
                    measured_idx_new(measured_idx_new == candidates(1)) = nj;
                end
                for kk = 2:(nb_merging + nj_is_measured)
                    measured_idx_new(measured_idx_new == candidates(kk)) = [];
                end

                % Condition to check that if we do the merge, we don't go
                % below the bound of nbsignals
                nbsignals_new = length(excited_idx) + length(measured_idx_new);
                if nbsignals_new >= nbsignals_bound

                    % Condition to check that if we do the merge, we still
                    % have nb_exc * nb_meas >= nb_unknowns
                    nb_equations = length(excited_idx) * length(measured_idx_new);
                    if nb_equations >= nbunknowns

                        % Checking that merging does not create meeting nodes in other subgraph.
                        % For each measurement in candidates, collect its assigned unknown edges.
                        % For each of those unknown edges, get the source node.
                        % Then compute rank T(source nodes, B).
                        assigned_unknowns = (sigma_C == candidates(1));
                        for mm = 2:(nb_merging + nj_is_measured)
                            % | = Binary OR
                            assigned_unknowns = assigned_unknowns | (sigma_C == candidates(mm));
                        end
                        sources = sources_unknowns(assigned_unknowns);
                        if rank(T_rand(sources, excited_idx)) == numel(sources)

                            % Construction of assigned_B :
                            % assigned_B collects all excitations assigned to all
                            % unknown edges assigned to each candidate measurement
                            nb_assigned_B = sum(ismember(sigma_C, candidates));
                            assigned_B = zeros(nb_assigned_B, 1);
                            idx = 1;
                            for mm = 1:(nb_merging + nj_is_measured)
                                edges_assigned = find(sigma_C == candidates(mm));
                                nb_edges_assigned = length(edges_assigned);
                                for nn = 1:nb_edges_assigned
                                    assigned_B(idx) = sigma_B(edges_assigned(nn));
                                    idx = idx + 1;
                                end
                            end

                            % Verifying that sigma remains bijective : all 
                            % corresponding unknown edges must have different excitation
                            % if numel(assigned_B) == numel(unique(assigned_B))
                            % Same but more efficient :
                            if all(diff(sort(assigned_B)) ~= 0)

                                % For each bi in assigned_B, check the vdp
                                % condition between the ending nodes of
                                % the unknown edges assigned to bi and C.
                                condition_vdp_bi = true;
                                for bi = assigned_B'
                                    assigned_unknowns_bi = (sigma_B == bi);
                                    targets_bi = targets_unknowns(assigned_unknowns_bi);
                                    rank_bi = rank(T_rand(measured_idx_new, targets_bi));
                                    if rank_bi < numel(targets_bi)
                                        condition_vdp_bi = false;
                                    end
                                end
                                if condition_vdp_bi

                                    % Identifiability test with rank(K)
                                    % Smart indexing for K
                                    nB = length(excited_idx);
                                    nC = length(measured_idx_new);
                                    indices = zeros(nB*nC,1);
                                    for l = 1:nB*nC
                                        indices(l) = (excited_idx(ceil(l/nC)) - 1)*L ...
                                            + measured_idx_new(mod(l-1,nC) + 1);
                                    end

                                    % If rank(K) = nbunknowns, do the merge :
                                    % Replace candidates by nj in sigma_C
                                    if rank(Kfull(indices,:)) == nbunknowns
                                        sigma_C(sigma_C == candidates(1)) = nj;
                                        for jj = 2:nb_merging
                                            sigma_C(sigma_C == candidates(jj)) = nj;
                                        end

                                        % Recursive call
                                        if ~B_done
                                            [excited_idx, measured_idx_out] = ...
                                                merging_uniform_B(graph, unknowns, ...
                                                targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
                                                excited_idx, measured_idx_new, sigma_B, sigma_C, stop, 0);
                                            return;
                                        else
                                            [excited_idx, measured_idx_out] = ...
                                                merging_uniform_C(graph, unknowns, ...
                                                targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
                                                excited_idx, measured_idx_new, sigma_B, sigma_C, stop, B_done);
                                            return;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
% We tried all possibilities of measurement merging => We set C_done to 1.
% If B_done == 0, we try merging excitations
if ~B_done
    [excited_idx, measured_idx_out] = merging_uniform_B(graph, unknowns, ...
        targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
        excited_idx, measured_idx, sigma_B, sigma_C, stop, 1);
else % Otherwise, return the current solution
    measured_idx_out = measured_idx;
end

end