function[Results_banks,Results,assetprices,FailCount,cum_Fails,capitalshortfall,num_ActiveBanks,num_ActiveAssets,...
    numnodes,numedges,density,avdegree,mu_A,mu_B,CB_TOTallotment,IBN_adjmat,OPN_adjmat,NNT_matrices,NMT_matrices,meanfail_index] = ...
    simcollect(abmresults,abmresults_label,n_sims,sim_Results_banks,sim_Results_agg,sim_Results_av,sim_Results_min,sim_Results_max,...
    sim_assetprices,sim_FailCount,sim_num_ActiveBanks,sim_num_ActiveAssets,...
    sim_numnodes,sim_numedges,sim_density,sim_avdegree,sim_mu_A,sim_mu_B,...
    sim_IBN_adjmat,sim_OPN_adjmat,sim_NNT_matrices,sim_NMT_matrices,T,sim_CB_TOTallotment)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average results across simulations before passing onto visualisation functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Averaging across ABM runs to obtain balance sheet, interbank market and firesale results for visualitation
Results(:,:,1) = nanmean(sim_Results_av,3);
Results(:,:,2) = nanmean(sim_Results_min,3);
Results(:,:,3) = nanmean(sim_Results_max,3);
Results(:,:,4) = nanmean(sim_Results_agg,3);

assetprices = nanmean(sim_assetprices,3);

% Average vector of bank failures across simulations
FailCount(1,:) = mean(sim_FailCount);
FailCount(2,:) = min(sim_FailCount);
FailCount(3,:) = max(sim_FailCount);

% Cumulative vector of failed banks after averaging across runs
cum_Fails(1,:) = cumsum(FailCount(1,:));
cum_Fails(2,:) = cumsum(FailCount(2,:));
cum_Fails(3,:) = cumsum(FailCount(3,:));

NumFails             = zeros(1,n_sims);
sim_capitalshortfall = zeros(n_sims,T);

% Identity of failed banks varies across ABM runs: Compute total capital loss due to failures in each run
for k = 1:n_sims
    Failures_sims        = abmresults.(abmresults_label{k}).FailedBanks; 
    NumFails(k)          = numel(Failures_sims); 
    capitalshortfall_vec = zeros(NumFails(k),T);

    for i=1:NumFails(k)
        capitalshortfall_vec(Failures_sims(i),abmresults.(abmresults_label{k}).vars(Failures_sims(i)).failtime) = ...
            abmresults.(abmresults_label{k}).vars(Failures_sims(i)).balancesheet{'Capital',end};
    end   
    sim_capitalshortfall(k,:) = abs(sum(capitalshortfall_vec));        
end

capitalshortfall(1,:) = mean(sim_capitalshortfall);
capitalshortfall(2,:) = min(sim_capitalshortfall);
capitalshortfall(3,:) = max(sim_capitalshortfall);

% Number of active banks and assets

num_ActiveBanks(1,:) = mean(sim_num_ActiveBanks);
num_ActiveBanks(2,:) = min(sim_num_ActiveBanks);  
num_ActiveBanks(3,:) = max(sim_num_ActiveBanks);  

num_ActiveAssets(1,:) = mean(sim_num_ActiveAssets);
num_ActiveAssets(2,:) = min(sim_num_ActiveAssets);  
num_ActiveAssets(3,:) = max(sim_num_ActiveAssets); 


% Network structure evolution across ABM runs
numnodes(1,:) = mean(sim_numnodes);
numnodes(2,:) = min(sim_numnodes);
numnodes(3,:) = max(sim_numnodes);

numedges(1,:) = mean(sim_numedges);
numedges(2,:) = min(sim_numedges);
numedges(3,:) = max(sim_numedges);

density(1,:)  = mean(sim_density);
density(2,:)  = min(sim_density);
density(3,:)  = max(sim_density);

avdegree(1,:) = mean(sim_avdegree);
avdegree(2,:) = min(sim_avdegree);
avdegree(3,:) = max(sim_avdegree);   

mu_A(1,:) = mean(sim_mu_A);
mu_A(2,:) = min(sim_mu_A);
mu_A(3,:) = max(sim_mu_A);

mu_B(1,:) = mean(sim_mu_B);
mu_B(2,:) = min(sim_mu_B);
mu_B(3,:) = max(sim_mu_B);

CB_TOTallotment(1,:) = mean(sim_CB_TOTallotment(:,:,1)); % Borrower refinancing
CB_TOTallotment(2,:) = mean(sim_CB_TOTallotment(:,:,2)); % Lender refinancing
CB_TOTallotment(3,:) = mean(sim_CB_TOTallotment(:,:,3)); % FRFA refinancing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamic adjacency matrices across simulations and policy scenarios
%% --> To be passed onto network visualisation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each policy experiment (S=1,2), choose ABM simulation in which total number of failures
% is closest to the average value across all simulations

mean(NumFails)
NumFails

[~,meanfail_index] = min(abs(NumFails-mean(NumFails)));

IBN_adjmat = sim_IBN_adjmat(:,:,:,meanfail_index);
OPN_adjmat = sim_OPN_adjmat(:,:,:,meanfail_index);

NNT_matrices = sim_NNT_matrices(:,:,:,:,meanfail_index);
NMT_matrices = sim_NMT_matrices(:,:,:,:,meanfail_index);

Results_banks = sim_Results_banks(:,:,:,meanfail_index);

end