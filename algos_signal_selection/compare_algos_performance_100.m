% Runs our algos on all networks with 100 nodes, for 10 sets
% of unknowns, 10 times each.
% Returns the average performances of each algo for each network.
% Prints those performances on the command line, and saves them in a file.

% Took 5h to run with
% ntrials = 10;
% nsamples = 10;
% npass_all = 10;

function performance = compare_algos_performance_100

start = tic;

filenames = ["ieee118", "watts100", "lattice10x10", "erdos100", "rgg100"];

npass_simug = 10;
npass_pseudotree = 10;
npass_merging = 10;
npass_pruning = 10;

% Number of sets of unknowns for each network
ntrials = 10;

% Number of samples for each set of unknowns of each network
% i.e. number of samples for each algo experiment 
nsamples = 10;

nbnets = length(filenames);
network = filenames';
avg_degree = zeros(nbnets,1);
stop_vec = zeros(nbnets,1);

nbsignals_simug = zeros(nbnets, 1);
nbsignals_pseudotree = zeros(nbnets, 1);
nbsignals_merging = zeros(nbnets, 1);
nbsignals_pruning = zeros(nbnets, 1);
nbsignals_exhaustive = zeros(nbnets, 1);
time_simug = zeros(nbnets, 1);
time_pseudotree = zeros(nbnets, 1);
time_merging = zeros(nbnets, 1);
time_pruning = zeros(nbnets, 1);
time_exhaustive = zeros(nbnets, 1);

% Loop on all networks
for j = 1:nbnets
    fullname = filenames(j) + ".mat";
    data = load(fullname);
    graph = data.graph;
    unknowns_tensor = data.unknowns_tensor;
    L = size(graph,1);
    nbedges = sum(graph, 'all');
    avg_degree(j) = 2 * nbedges / L;

    [nrun, stop] = parameters_merging(data.topology, L);
    stop_vec(j) = stop;

    nbsignals_simug_list = zeros(ntrials,1);
    nbsignals_pseudotree_list = zeros(ntrials,1);
    nbsignals_merging_list = zeros(ntrials,1);
    nbsignals_pruning_list = zeros(ntrials,1);
    nbsignals_exhaustive_list = zeros(ntrials,1);
    time_simug_list = zeros(ntrials,1);
    time_pseudotree_list = zeros(ntrials,1);
    time_merging_list = zeros(ntrials,1);
    time_pruning_list = zeros(ntrials,1);
    time_exhaustive_list = zeros(ntrials,1);

    % Loop on all set of unknowns for network j
    for i = 1:ntrials

        nbsignals_simug_samples = zeros(nsamples,1);
        nbsignals_pseudotree_samples = zeros(nsamples,1);
        nbsignals_merging_samples = zeros(nsamples,1);
        nbsignals_pruning_samples = zeros(nsamples,1);
        nbsignals_exhaustive_samples = zeros(nsamples,1);
        time_simug_samples = zeros(nsamples,1);
        time_pseudotree_samples = zeros(nsamples,1);
        time_merging_samples = zeros(nsamples,1);
        time_pruning_samples = zeros(nsamples,1);
        time_exhaustive_samples = zeros(nsamples,1);

        % Each algo is run nsamples times
        for k = 1:nsamples
            tic
            [excited_simug, measured_simug] = simug(graph, unknowns_tensor(:,:,i), npass_simug);
            time_simug_samples(k) = toc;
            nB_simug = nnz(excited_simug);
            nC_simug = nnz(measured_simug);
            nbsignals_simug_samples(k) = nB_simug + nC_simug;
            
            tic
            [excited_pseudotree, measured_pseudotree] = pseudotree(graph, unknowns_tensor(:,:,i), npass_pseudotree);
            time_pseudotree_samples(k) = toc;
            nB_pseudotree = nnz(excited_pseudotree);
            nC_pseudotree = nnz(measured_pseudotree);
            nbsignals_pseudotree_samples(k) = nB_pseudotree + nC_pseudotree;

            tic
            [excited_merging, measured_merging] = merging(graph, unknowns_tensor(:,:,i), nrun, npass_merging, stop);
            time_merging_samples(k) = toc;
            nB_merging = nnz(excited_merging);
            nC_merging = nnz(measured_merging);
            nbsignals_merging_samples(k) = nB_merging + nC_merging;

            tic
            [excited_pruning, measured_pruning] = pruning(graph, unknowns_tensor(:,:,i), npass_pruning);
            time_pruning_samples(k) = toc;
            nB_pruning = nnz(excited_pruning);
            nC_pruning = nnz(measured_pruning);
            nbsignals_pruning_samples(k) = nB_pruning + nC_pruning;

            if L <= 14
                tic
                [excited_exhaustive, measured_exhaustive] = exhaustive(graph, unknowns_tensor(:,:,i));
                time_exhaustive_samples(k) = toc;
                nB_exhaustive = nnz(excited_exhaustive);
                nC_exhaustive = nnz(measured_exhaustive);
                nbsignals_exhaustive_samples(k) = nB_exhaustive + nC_exhaustive;
            end
        end

        nbsignals_simug_list(i) = mean(nbsignals_simug_samples);
        nbsignals_pseudotree_list(i) = mean(nbsignals_pseudotree_samples);
        nbsignals_merging_list(i) = mean(nbsignals_merging_samples);
        nbsignals_pruning_list(i) = mean(nbsignals_pruning_samples);
        nbsignals_exhaustive_list(i) = mean(nbsignals_exhaustive_samples);
        time_simug_list(i) = mean(time_simug_samples);
        time_pseudotree_list(i) = mean(time_pseudotree_samples);
        time_merging_list(i) = mean(time_merging_samples);
        time_pruning_list(i) = mean(time_pruning_samples);
        time_exhaustive_list(i) = mean(time_exhaustive_samples);

    end

    nbsig_sim = mean(nbsignals_simug_list);
    nbsig_ptree = mean(nbsignals_pseudotree_list);
    nbsig_merg = mean(nbsignals_merging_list);
    nbsig_prun = mean(nbsignals_pruning_list);
    nbsig_exhau = mean(nbsignals_exhaustive_list);
    time_sim = mean(time_simug_list);
    time_ptree = mean(time_pseudotree_list);
    time_merg = mean(time_merging_list);
    time_prun = mean(time_pruning_list);
    time_exhau = mean(time_exhaustive_list);

    nbsignals_simug(j) = nbsig_sim;
    nbsignals_pseudotree(j) = nbsig_ptree;
    nbsignals_merging(j) = nbsig_merg;
    nbsignals_pruning(j) = nbsig_prun;
    nbsignals_exhaustive(j) = nbsig_exhau;
    time_simug(j) = time_sim;
    time_pseudotree(j) = time_ptree;
    time_merging(j) = time_merg;
    time_pruning(j) = time_prun;
    time_exhaustive(j) = time_exhau;

    disp(filenames(j));
    performance_j = table(nbsig_sim, nbsig_ptree, nbsig_merg, nbsig_prun, ...
        nbsig_exhau, time_sim, time_ptree, time_merg, time_prun, time_exhau);
    disp(performance_j);

    performance = table(network, avg_degree, nbsignals_simug, nbsignals_pseudotree, nbsignals_merging, ...
        nbsignals_pruning, nbsignals_exhaustive, ...
    time_simug, time_pseudotree, time_merging, time_pruning, time_exhaustive, ...
    stop_vec);
    t = char(datetime('now','Format','MM-dd''-T:''HH:mm'));
    writetable(performance,['performance_sofar_100_', t, '.txt']);
end

performance = table(network, avg_degree, nbsignals_simug, nbsignals_pseudotree, nbsignals_merging, ...
    nbsignals_pruning, nbsignals_exhaustive, ...
    time_simug, time_pseudotree, time_merging, time_pruning, time_exhaustive, ...
    stop_vec);
t = char(datetime('now','Format','MM-dd''-T:''HH:mm'));
writetable(performance,['performance_all_100_', t, '.txt']);

elapsed = toc(start)
end