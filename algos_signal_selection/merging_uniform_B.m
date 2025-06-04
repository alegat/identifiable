% Merging Excitations - Uniform pass

function [excited_idx_out, measured_idx] = merging_uniform_B(graph, unknowns, ...
    targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
    excited_idx, measured_idx, sigma_B, sigma_C, stop, C_done)

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
    excited_idx_out = excited_idx;
    return
end

% Loop on all excitations, in random order
excited_idx = excited_idx(randperm(length(excited_idx)));
for i = 1:length(excited_idx)
    bi = excited_idx(i);
    inneighbours = find(ancestors(bi,:));
    inneighbours = inneighbours(randperm(length(inneighbours)));

    % Loop on all multi-generational inneighbours of excitation bi
    for j = 1:length(inneighbours)
        nj = inneighbours(j);
        nj_is_excited = ismember(nj, excited_idx);
        outneighbours = find(ancestors(:,nj));
        outneighbours = outneighbours(randperm(length(outneighbours)));
        intersection = intersect(outneighbours, excited_idx);

        % To avoid trying to merge too many excitations that would go below
        % the bound of nbsignals
        nb_merging_max = length(excited_idx) + length(measured_idx) ...
            - nbsignals_bound + 1 - nj_is_excited;
        nb_merging_start = min(nb_merging_max, length(intersection));

        % nb_merging is the nb of excitations we try to merge into nj
        % We are greedy: we start with the highest possible nb_merging and
        % then try with smaller ones
        for nb_merging = nb_merging_start:-1:(2 - nj_is_excited)

            if nb_merging == 1
                groups = intersection;
            else
                groups = nchoosek(intersection,nb_merging);
            end
            nb_groups = size(groups, 1);
            groups = groups(randperm(nb_groups),:);

            % Loop on all groups of excitations computed by nchoosek
            for k = 1:nb_groups
                candidates = groups(k,:);
                if nj_is_excited
                    candidates = [nj, candidates];
                end

                % Computation of excited_idx_new
                excited_idx_new = excited_idx;
                if ~ nj_is_excited
                    excited_idx_new(excited_idx_new == candidates(1)) = nj;
                end
                for kk = 2:(nb_merging + nj_is_excited)
                    excited_idx_new(excited_idx_new == candidates(kk)) = [];
                end

                % Condition to check that if we do the merge, we don't go
                % below the bound of nbsignals
                nbsignals_new = length(excited_idx_new) + length(measured_idx);
                if nbsignals_new >= nbsignals_bound

                    % Condition to check that if we do the merge, we still
                    % have nb_exc * nb_meas >= nb_unknowns
                    nb_equations = length(excited_idx_new) * length(measured_idx);
                    if nb_equations >= nbunknowns

                        % Checking that merging does not create meeting nodes in other subgraph.
                        % For each excitation in candidates, collect its assigned unknown edges.
                        % For each of those unknown edges, get the target node.
                        % Then compute rank T(C, target nodes).
                        assigned_unknowns = (sigma_B == candidates(1));
                        for ll = 2:(nb_merging + nj_is_excited)
                            % | = Binary OR
                            assigned_unknowns = assigned_unknowns | (sigma_B == candidates(ll));
                        end
                        targets = targets_unknowns(assigned_unknowns);
                        if rank(T_rand(measured_idx, targets)) == numel(targets)

                            % assigned_C collects all measurements assigned to all
                            % unknown edges assigned to each candidate excitation
                            nb_assigned_C = sum(ismember(sigma_B, candidates));
                            assigned_C = zeros(nb_assigned_C, 1);
                            idx = 1;
                            for mm = 1:(nb_merging + nj_is_excited)
                                edges_assigned = find(sigma_B == candidates(mm));
                                nb_edges_assigned = length(edges_assigned);
                                for nn = 1:nb_edges_assigned
                                    assigned_C(idx) = sigma_C(edges_assigned(nn));
                                    idx = idx + 1;
                                end
                            end

                            % Verifying that sigma remains bijective = all
                            % corresponding unknown edges must have different measurement
                            % if numel(assigned_C) == numel(unique(assigned_C))
                            % Same but more efficient :
                            if all(diff(sort(assigned_C)) ~= 0)

                                % For each ci in assigned_C, check the vdp
                                % condition between the starting nodes of
                                % the unknown edges assigned to ci and B.
                                condition_vdp_ci = true;
                                for ci = assigned_C'
                                    assigned_unknowns_ci = (sigma_C == ci);
                                    sources_ci = sources_unknowns(assigned_unknowns_ci);
                                    rank_ci = rank(T_rand(sources_ci, excited_idx_new));
                                    if rank_ci < numel(sources_ci)
                                        condition_vdp_ci = false;
                                    end
                                end
                                if condition_vdp_ci

                                    % Identifiability test with rank(K)
                                    % Smart indexing for K
                                    nB = length(excited_idx_new);
                                    nC = length(measured_idx);
                                    indices = zeros(nB*nC,1);
                                    for l = 1:nB*nC
                                        indices(l) = (excited_idx_new(ceil(l/nC)) - 1)*L ...
                                            + measured_idx(mod(l-1,nC) + 1);
                                    end

                                    % If rank(K) = nbunknowns, do the merge :
                                    % Replace candidates by nj in sigma_B
                                    if rank(Kfull(indices,:)) == nbunknowns
                                        sigma_B(sigma_B == candidates(1)) = nj;
                                        for jj = 2:nb_merging
                                            sigma_B(sigma_B == candidates(jj)) = nj;
                                        end

                                        % Recursive call
                                        if ~C_done
                                            [excited_idx_out, measured_idx] = ...
                                                merging_uniform_C(graph, unknowns, ...
                                                targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
                                                excited_idx_new, measured_idx, sigma_B, sigma_C, stop, 0);
                                            return;
                                        else
                                            [excited_idx_out, measured_idx] = ...
                                                merging_uniform_B(graph, unknowns, ...
                                                targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
                                                excited_idx_new, measured_idx, sigma_B, sigma_C, stop, C_done);
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
% We tried all possibilities of excitation merging => We set B_done to 1.
% If C_done == 0, we try merging measurements
if ~C_done
    [excited_idx_out, measured_idx] = merging_uniform_C(graph, unknowns, ...
        targets_unknowns, sources_unknowns, ancestors, T_rand, Kfull, ...
        excited_idx, measured_idx, sigma_B, sigma_C, stop, 1);
else % Otherwise, return the current solution
    excited_idx_out = excited_idx;
end

end