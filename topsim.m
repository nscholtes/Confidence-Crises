function [Networkdata,Graphdata,ABMresults] = topsim()

%--------------------------------------------------------------------------
% MATLAB code for 'Exploring the counterparty-liquidity risk nexus using agent-based network model of the interbank market
% by Nicolas K. Scholtes (2017)

% Top-level function containing all calibrated parameters. Calls the following lower-level functions:
%--------------------------------------------------------------------------
% NETWORK GENERATION
%--------------------------------------------------------------------------
%%% - netgen.m: Draw node fitness from a truncated power law
%%% - opnet.m: Generation of the random bipartite graph representing the
%%%   network of overlapping portfolios between banks
%--------------------------------------------------------------------------
% AGENT-BASED MODEL
%--------------------------------------------------------------------------
%%% 1. MAIN FUNCTIONS FOR ABM
%%%% - Phase1.m: Encodes  phase 1 - Exogenous deposit shock and  move to interbank market
%%%% - Phase2.m: Encodes phase 2 - Exogenous external asset shock, bank firesales and loan repayment
%%%% - Phase3.m: Encodes phase 3 - Bank insolvency computation and central bank intervention
%%%% - Periodsummary.m: Collects, summarises an organises results from each iteration step
%%% 2. INITIALIZING AND COLLECTING RESULTS
%%%% - siminitialize.m:
%%%% - simcollect.m:
%%% 3. RESULTS STORAGE AND VISUALISATION
%%%% - Results_balancesheet.m:
%%%% - Results_IBM.m:
%%%% - Results_securitiesmarket.m:
%%%% - Results_failures.m:
%%%% - abmresults.m:
%%%% - viewnetworks.m:

% Copyright (C) Nicolas K. Scholtes, 2017
% Distributed under GPL 3.0
%--------------------------------------------------------------------------
%% Script controls
%--------------------------------------------------------------------------
close all; clear; clc;

% Output figures to directory with LaTeX draft for automatic updating
fig_output = '/Users/nscholte/Desktop/Research/Ch.1 - Confidence crises/Drafts/Current version/Figures/';

% timestamp recording date of current simulation
timestamp = datestr(datetime('today'));

% Create directory for simulation logs (with timestamp)
if ~exist('Logs','dir') && ~exist('Logs/Network Generation','dir')
    mkdir Logs
    mkdir Logs 'Network Generation' 
end

fileID_D   = fopen(strcat('Logs/detailed_log_',timestamp,'.txt'),'w');
fileID_S   = fopen(strcat('Logs/summary_log_',timestamp,'.txt'),'w');
fileID_IBN = fopen(strcat('Logs/Network Generation/IBN_',timestamp,'.txt'),'w');
%fileID_OPN = fopen(strcat('Logs/Network Generation/OPN_',datestr(datetime('today')),'.txt'),'w');

% Select whether to create a detailed simulation log  (Input: Y/N)
writeoption = 'N';

tol = 1.e-6; % Tolerance value for == conditional statements

%--------------------------------------------------------------------------
% CALIBRATION OF SIMULATION PARAMETERS
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bank size/fitness distribution parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_banks     = 50;        % Number of banks
a_min       = 5;         % Minimum bank size
a_max       = 100;       % Maximum bank size
gamma       = 2;         % Power law exponent for bank size distribution
d           = 0.5;       % Calibrate network density

Pars_netgen = [d,a_min,a_max,gamma,n_banks];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlapping portfolio parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_assets = 50; % Number of external assets
av_div   = 10;   % Average number of assets held per bank (diversification)

Pars_opnet = [m_assets,av_div];

%--------------------------------------------------------------------------
% Agent-Based Model Parameters
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation and Model timing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_networks = 100;   % Number of simulated networks
n_sims     = 10;     % Number of ABM simulations

T_sim = 250;        % Number of iteration steps
T     = T_sim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial balance sheet weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.9;   % external asset/asset ratio
beta  = 0.92;  % deposit/liability ratio

Pars_balancesheet= [alpha beta];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liquidity and portfolio shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Liquidity (deposit shock)
theta = 0.5;   % Multiplier for mean-reverting component
SPF   = 0.025; % Shock proportionality factor

% Portfolio shock
shocktime   = T/2:(T/2)+50; 

num_shocked_ea_S = m_assets/10;
eashock_LB_S     = 0.9;
eashock_UB_S     = 1.1;
market_depth_S   = zeros(1,m_assets); % Market depth of assets used in firesale computation
    
num_shocked_ea_C = m_assets/2;
eashock_LB_C     = 0.75;
eashock_UB_C     = 1;
market_depth_C   = 0.2*ones(1,m_assets);
    
%Collecting parameters
Pars_dshock  = [theta SPF];
Pars_pshock = [num_shocked_ea_S, eashock_LB_S, eashock_UB_S market_depth_S;
    num_shocked_ea_C, eashock_LB_C, eashock_UB_C market_depth_C];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Policy experiment parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose simulation type: stable or crisis:
%%% Baseline: Small liquidity shocks and high market liquidity
%%% Crisis: Stable state followed by crisis (large shocks and low market liquidity) lasting a fixed amount of periods

MRR                  = 0.02;     % Minimum reserve requirement
refinancingfrequency = 5;        % Central bank refinancing frequency (weekly)
simtype              = 'crisis';
FRFAtime             = T/2+25:T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interest rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_z  = 0.05; % Illiquid investment rate
r_d  = 0.02; % Deposit rate
r_e  = 0.03; % External asset return
r_b  = 0.04; % Initial interbank rate

Pars_interestrates = [r_z r_d r_e r_b];

%--------------------------------------------------------------------------
% Generating the network
%--------------------------------------------------------------------------

a_sim          = zeros(n_banks,n_networks);
ibn_adjmat_sim = zeros(n_banks,n_banks,n_networks);
density_vec    = zeros(1,n_networks);

for k = 1:n_networks
    fprintf(1,'Network iteratiom: %d/%d\n',k,n_networks);
    [a_sim(:,k),ibn_adjmat_sim(:,:,k)] = netgen(Pars_netgen,fileID_IBN,fig_output);   
    density_vec(k) = density_und(ibn_adjmat_sim(:,:,k)); % Compute density of each simulated network
    clc
end

median_density = median(density_vec); % Compute median density

fprintf('The median density across %d simulated networks is %.3f\n',n_networks,median_density)

% Find index of median density
if ismember(median_density,density_vec)
    median_index = find(median_density == density_vec);
    if numel(median_index)>1
        median_index = datasample(median_index,1);
    end
else
    %median_index =  (min(find(median_density>median(density_vec))) + max(find(median_density<median(density_vec))))/2;
    [~,ord] = sort(density_vec);
    median_index = ord(floor(n_networks/2)+(rem(n_networks,2)==0):floor(n_networks/2)+1);
    %median_density = mean(density_vec(median_index));

end

% Use index to output the final interbank exposure adjacency matrix passed on to the ABM
a = a_sim(:,median_index);

fig_output_IBN   = strcat(fig_output,'Network/');

figure
hist(a);
title('Empirical distribution of node fitness values')
xlabel('a')
ylabel('f(a)')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBN,'sizedist.pdf'));

ibn_adjmat_init = ibn_adjmat_sim(:,:,median_index);

[opn_adjmat_init]   = opnet(n_banks,Pars_opnet,T,fig_output);

clearvars a_sim ibn_adjmat_sim density_vec  median_index ord

save('network.mat'); % Save workspace as .mat file

close all;
clc

%--------------------------------------------------------------------------
%% Agent-Based Model
%--------------------------------------------------------------------------

load('network.mat')

sim_vector = {'off','on'};

FIN_RESULTS          = zeros(15,T_sim,4,2);
FIN_ASSETPRICES      = zeros(4*T_sim,m_assets,2);
FIN_FAILCOUNT        = zeros(3,T_sim,2);
FIN_CUM_FAILS        = zeros(3,T_sim,2);
FIN_CAPITALSHORTFALL = zeros(3,T_sim,2);
FIN_NUM_ACTIVEBANKS  = zeros(3,T_sim,2);
FIN_NUM_ACTIVEASSETS = zeros(3,T_sim,2);
FIN_NUMNODES         = zeros(3,T_sim,2);
FIN_NUMEDGES         = zeros(3,T_sim,2);
FIN_DENSITY          = zeros(3,T_sim,2);
FIN_AVDEGREE         = zeros(3,T_sim,2);
FIN_MU_A             = zeros(3,T_sim,2);
FIN_MU_B             = zeros(3,T_sim,2);
FIN_CB_TOTALLOTMENT  = zeros(3,T_sim,2);

FIN_NNT_matrices     = zeros(n_banks,n_banks,T,8,2);
FIN_NMT_matrices     = zeros(n_banks,m_assets,4*T,2,2);

FIN_IBN_ADJMAT       = zeros(n_banks,n_banks,T_sim,2);
FIN_OPN_ADJMAT       = zeros(n_banks,m_assets,T_sim,2);

FIN_RESULTS_BANKS = zeros(n_banks,T_sim,14,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outer loop: Cycle through different policy and environment scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for S = 1:2

if strcmp(simtype,'baseline')
    CBintervention = sim_vector{S};
    FRFA           = 'NA';
    sim_label      = {'noCB','CB_noFRFA'};
elseif strcmp(simtype,'crisis')
    CBintervention = 'on';
    sim_label      = {'CB_noFRFA','CB_FRFA'};
end

T = T_sim;

sim_IBN_adjmat = zeros(n_banks,n_banks,T,n_sims);
sim_OPN_adjmat = zeros(n_banks,m_assets,T,n_sims);

sim_NNT_matrices = zeros(n_banks,n_banks,T,8,n_sims);
sim_NMT_matrices = zeros(n_banks,m_assets,4*T,2,n_sims);
sim_TM_matrices  = zeros(T,m_assets,2);

% Storing summarised matrix for each run of the ABM
sim_Results_banks = zeros(n_banks,T,14,n_sims);
sim_Results_agg   = zeros(15,T,n_sims);
sim_Results_av    = zeros(15,T,n_sims);
sim_Results_min   = zeros(15,T,n_sims);
sim_Results_max   = zeros(15,T,n_sims);

sim_assetprices   = ones(4*T,m_assets,n_sims);
sim_failcount     = zeros(n_sims,T);

sim_num_ActiveBanks  = zeros(n_sims,T);
sim_num_ActiveAssets = zeros(n_sims,T);

sim_numnodes  = zeros(n_sims,T);
sim_numedges  = zeros(n_sims,T);
sim_density   = zeros(n_sims,T);
sim_avdegree  = zeros(n_sims,T);
sim_mu_A      = zeros(n_sims,T);
sim_mu_B      = zeros(n_sims,T);

sim_CB_TOTallotment = zeros(n_sims,T,3);

close all 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested loop level 1: For each policy/environment scenario
% run the ABM $n_sims$ times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for k = 1:n_sims

% Define label to store cleaned ABM results for each run
abmresults_label{k} = strcat('sim',num2str(k));

% Function to initialize all vectors, matrices and structures for current simulation run
[ActiveBanks,ActiveAssets,banks,banksfail,num_ActiveBanks,num_FailedBanks,...
    num_ActiveAssets,num_InactiveAssets,FailCount_vec,FailedBankID_mat,FailedBanks,...
    numedges,mu_A,mu_B,assetprices,NT_matrices,NNT_matrices,NMT_matrices,TM_matrices,...
    IBN_adjmat,OPN_adjmat,Results_banks,Results_agg,Results_av,Results_min,Results_max,Results_dAgg,...
    total_firesales_vec,d_assetprices,d_firesales,CB_TOTallotment] =...
    siminitialize(T_sim,n_banks,m_assets,ibn_adjmat_init,opn_adjmat_init);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested loop level 2: Run the ABM over the predetermiend number of iteration steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for t=1:T
        if t==1
            x = 1;
        else  
            x = t-1;
        end
        if strcmp(simtype,'crisis')
            if ~ismember(t,shocktime)
                Pars_pshock_current = Pars_pshock(1,:);
            elseif ismember(t,shocktime)
                Pars_pshock_current = Pars_pshock(2,:);
            end  
            if S == 1    % 1st crisis experiment: central bank intervenes normally but no FRFA
                FRFA = 'off';
            elseif S == 2 % 2nd crisis experiment: normal interventions and turn on FRFA at a predetermined time
                if ~ismember(t,FRFAtime)
                    FRFA = 'off';
                elseif ismember(t,FRFAtime)
                    FRFA = 'on';
                end
            end
        elseif strcmp(simtype,'baseline')
            Pars_pshock_current = Pars_pshock(1,:);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Print information on the current simulation type
        disp('--------------- SIMULATION PARAMETERS ---------------')
        fprintf('The current simulation type is: %s\n',simtype)
        if strcmp(simtype,'crisis')
            fprintf('--> Crisis starts at period %d and lasts %d periods\n',shocktime(1),shocktime(end)-shocktime(1))
        end
        fprintf('Central bank intervention is %s\n',CBintervention)
        fprintf('--> Fixed rate full allotment policy is %s\n',FRFA)
        if strcmp(FRFA,'on')
            disp('--> FRFA policy intervention active')
        end
        disp('-----------------------------------------------------')
        fprintf(1,'Current simulation: %d/%d\n',k,n_sims);
        fprintf(1,'--> Current period: %d/%d\n',t,T_sim);
% Phase 0: Initialisation, Population of the network and carrying over of variables from previous iterations
        [banks,NMT_matrices,assetprices] =...
            Phase0(banks,ActiveBanks,a,Pars_balancesheet,T,t,IBN_adjmat(:,:,x),NMT_matrices,OPN_adjmat(:,:,x),assetprices);
% Phase 1: Deposit shock, investment decision, interbank market
        [banks,DBV,DLV,NT_matrices,NNT_matrices] = ...
            Phase1(banks,ActiveBanks,IBN_adjmat(:,:,x),n_banks,beta,Pars_interestrates,Pars_dshock,MRR,...
            NT_matrices,NNT_matrices,t,fileID_D,writeoption,tol);
% Phase 2: Asset price shock, loan repayment, firesales and second round effects
        [banks,assetprices,OPN_adjmat(:,:,t),NMT_matrices,NNT_matrices,TM_matrices] = ...
            Phase2(banks,n_banks,ActiveBanks,ActiveAssets,OPN_adjmat(:,:,x),Pars_pshock_current,MRR,...
            NMT_matrices,NNT_matrices,TM_matrices,assetprices,DBV,DLV,t,fileID_D,writeoption,tol);
% Phase 3: Removal of insolvent banks and (possible) central bank intervention
        [banks,NMT_matrices,banksfail,n_banks,IBN_adjmat(:,:,t),OPN_adjmat(:,:,t),...
            FailCount_vec(t),FailedBankID_mat(t,:),FailedBanks,ActiveBanks,num_ActiveBanks(t),...
            num_FailedBanks(t),ActiveAssets,num_ActiveAssets(t),num_InactiveAssets(t),CB_TOTallotment(:,t)] =  ...
            Phase3(banks,DBV,DLV,NMT_matrices,banksfail,n_banks,m_assets,Pars_interestrates(1),...
            IBN_adjmat(:,:,x),OPN_adjmat(:,:,x),assetprices(4*t,:),...
            FailCount_vec(t),FailedBankID_mat(t,:),FailedBanks,ActiveBanks,ActiveAssets,...
            CBintervention,refinancingfrequency,FRFA,CB_TOTallotment(:,t),t,T,tol); 
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Current iteration summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        [Results_banks(:,t,:),Results_agg(:,1:t),Results_av(:,1:t),Results_min(:,1:t),Results_max(:,1:t),...
            Results_dAgg(:,t),total_firesales_vec(:,:,t),d_assetprices((2*t)-1:2*t,:),d_firesales(t),numedges(t),mu_A(t),mu_B(t)] =...
            Periodsummary(banks,ActiveBanks,num_ActiveBanks(t),num_ActiveAssets(t),n_banks,NNT_matrices(:,:,t,8),Results_banks(:,t,:),...
            Results_agg(:,1:t),Results_av(:,1:t),Results_min(:,1:t),Results_max(:,1:t),Results_dAgg(:,t),...
            m_assets,assetprices,total_firesales_vec(:,:,t),d_assetprices((2*t)-1:2*t,:),d_firesales(t),...
            IBN_adjmat(:,:,t),OPN_adjmat(:,:,t),FailCount_vec(t),FailedBankID_mat(t,:),DBV,DLV,t,fileID_S);  
        clc;
        if numel(ActiveBanks) == 0
            fprintf(1,'All banks failed by period %d\n',t)
            T = t;
            break
        end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intra-period housekeeping and collecting results across simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    [banks] = housekeeping(banks,n_banks);
    clearvars DBV DLV
    end
    
    % Adjacency matrix dynamics. DIM = N x N x T x K x 2. DIM(4) = number of ABM simulations   
    sim_NNT_matrices(:,:,:,:,k) = NNT_matrices;
    sim_NMT_matrices(:,:,:,:,k) = NMT_matrices;
    sim_TM_matrices(:,:,:,k)    = TM_matrices;
    
    sim_IBN_adjmat(:,:,:,k) = IBN_adjmat;
    sim_OPN_adjmat(:,:,:,k) = OPN_adjmat;
    
    % Balance sheet and interbank market dynamics. DIM = 15 x T x K
    sim_Results_banks(:,:,:,k)  = Results_banks;
    sim_Results_agg(:,:,k)      = Results_agg;
    sim_Results_av(:,:,k)       = Results_av;
    sim_Results_min(:,:,k)      = Results_min;
    sim_Results_max(:,:,k)      = Results_max;
    sim_assetprices(:,:,k)      = assetprices;
        
    % Network structure dynamics
    sim_failcount(k,:) = FailCount_vec;
    sim_numnodes(k,:)  = n_banks*ones(1,T_sim) - cumsum(FailCount_vec);
    sim_numedges(k,:)  = numedges;
    
    sim_num_ActiveBanks(k,:)  = num_ActiveBanks;
    sim_num_ActiveAssets(k,:) = num_ActiveAssets;
    
    sim_density(k,:)  = sim_numedges(k,:)./(sim_numnodes(k,:).*(sim_numnodes(k,:)-ones(1,T_sim)));
    sim_avdegree(k,:) = sim_numedges(k,:)./sim_numnodes(k,:);
    
    sim_mu_A(k,:) = mu_A;
    sim_mu_B(k,:) = mu_B;
    
    sim_CB_TOTallotment(k,:,1) = CB_TOTallotment(1,:);
    sim_CB_TOTallotment(k,:,2) = CB_TOTallotment(2,:);
    sim_CB_TOTallotment(k,:,3) = CB_TOTallotment(3,:);
    
    FIN_abmresults.(sim_label{S}).(abmresults_label{k}).vars         = abmresults(banks,n_banks,MRR);
    FIN_abmresults.(sim_label{S}).(abmresults_label{k}).ActiveBanks  = ActiveBanks;
    FIN_abmresults.(sim_label{S}).(abmresults_label{k}).FailedBanks  = FailedBanks; 
    FIN_abmresults.(sim_label{S}).(abmresults_label{k}).ActiveAssets = ActiveAssets;
    
clearvars banks banksfail abmresults Results_banks Results_agg Results_av Results_min Results_max Results_d_Agg...
    assetprices FailCount_vec FailedBankID_mat numedges ActiveBanks num_ActiveBanks FailedBanks num_FailedBanks...
    ActiveAssets num_ActiveAssets num_InactiveAssets mu_A mu_B CB_TOTallotment total_firesales_vec...
    IBN_adjmat OPN_adjmat NNT_matrices NMT_matrices NT_matrices TM_matrices

clc
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COLLECTING SIMULUATION RESULTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[FIN_RESULTS_BANKS(:,:,:,S),FIN_RESULTS(:,:,:,S) ,FIN_ASSETPRICES(:,:,S),FIN_FAILCOUNT(:,:,S),FIN_CUM_FAILS(:,:,S),...
    FIN_CAPITALSHORTFALL(:,:,S),FIN_NUM_ACTIVEBANKS(:,:,S),FIN_NUM_ACTIVEASSETS(:,:,S),...
    FIN_NUMNODES(:,:,S),FIN_NUMEDGES(:,:,S)  ,FIN_DENSITY(:,:,S),FIN_AVDEGREE(:,:,S),...
    FIN_MU_A(:,:,S),FIN_MU_B(:,:,S) ,FIN_CB_TOTALLOTMENT(:,:,S),...
    FIN_IBN_ADJMAT(:,:,:,S),FIN_OPN_ADJMAT(:,:,:,S),FIN_NNT_matrices(:,:,:,:,S),FIN_NMT_matrices(:,:,:,:,S),meanfail_index(S)] =...
    simcollect(FIN_abmresults.(sim_label{S}),abmresults_label,n_sims,...
    sim_Results_banks,sim_Results_agg,sim_Results_av,sim_Results_min,sim_Results_max,...
    sim_assetprices,sim_failcount,sim_num_ActiveBanks,sim_num_ActiveAssets,sim_numnodes,sim_numedges,...
    sim_density,sim_avdegree,sim_mu_A,sim_mu_B,sim_IBN_adjmat,sim_OPN_adjmat,sim_NNT_matrices,sim_NMT_matrices,...
    T_sim,sim_CB_TOTallotment);

clearvars sim_Results_banks sim_Results_agg sim_Results_av sim_Results_min sim_Results_max sim_assetprices sim_failcount...   
sim_num_ActiveBanks sim_num_ActiveAssets sim_numnodes sim_numedges sim_density sim_avdegree sim_mu_A...     
sim_mu_B sim_CB_TOTallotment 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLEAN UP  AND SAVE WORKSPACE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
clearvars a a_max a_min alpha av_div beta d_assetprices d_firesales eashock_LB_C eashock_LB_S...
    eashock_UB_C eashock_UB_S fileID_D fileID_IBN fileID_S gamma dTOT_matrices market_depth_C...
    market_depth_S num_shocked_ea_C num_shocked_ea_S Pars_balancesheet Pars_dshock Pars_interestrates...
    Pars_netgen Pars_opnet Pars_policy Pars_pshock Pars_pshock_current psi r_b r_d r_e r_z SPF theta ...
    
save(char(strcat('output_',simtype,'_',sim_label(S),'_',timestamp,'.mat'))); % Save workspace as .mat file

%--------------------------------------------------------------------------
%% Model output
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NETWORK VISUALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fix parameters used for network visualization
timeslices_IBN = [120 140 160 180];
timeslices_OPN = [160 165 170 175 180 200];
NNTvisualize   = [1 5 8];  % Select which weighted matrices to visualize: (1 = ex-post adjacency matrix, 5 = Loan matrix, 8 = Interbank rate matrix)
NTvisualize    = 1;        % Total assets of active banks in each period
markernorm     = 10;       % Normalization factor such that marker size of largest node = $markernorm$
edgewnorm      = 10;       % Normalization factor such that width of largest edge = %edgewnorm$
graphlayout    = 'circle'; % Choose from:  'auto' (default) | 'circle' | 'force' | 'layered' | 'subspace' | 'force3' | 'subspace3'

networkparameters = [markernorm edgewnorm graphlayout];

[Networkdata,Graphdata] = viewnetworks(fig_output,networkparameters,FIN_RESULTS_BANKS(:,[1 timeslices_IBN],NTvisualize,1),...
    [1 timeslices_IBN],[1 timeslices_OPN],...
    cat(3,ibn_adjmat_init,FIN_IBN_ADJMAT(:,:,timeslices_IBN,1)),...
    FIN_NNT_matrices(:,:,[1 timeslices_IBN],NNTvisualize,1),...
    cat(3,opn_adjmat_init,FIN_OPN_ADJMAT(:,:,timeslices_OPN,1)),tol);
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL OUTPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalizegraph = 'NormOff_';  % NormOn_/NormOff_ to normalize balance sheet variables by total assets
results_format = 'av_';       % Choose from: agg_ | av_. Determines how results are presented
includebounds  = 'BOff';     % BOn_/BOff_: Choose whether to include upper and lower bound results across ABM runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BALANCE SHEET DYNAMICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BS_outputlabel = strcat('[',normalizegraph,results_format,includebounds,']');
Results_balancesheet(fig_output,simtype,FIN_RESULTS,...
    T_sim,shocktime,FRFAtime,normalizegraph,results_format,includebounds,BS_outputlabel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERBANK MARKET DYNAMICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IBM_outputlabel = strcat('[',results_format,includebounds,']');
Results_IBM(fig_output,simtype,FIN_RESULTS,FIN_CB_TOTALLOTMENT,T_sim,shocktime,FRFAtime,...
    results_format,includebounds,IBM_outputlabel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECURITIES MARKET DYNAMICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SM_outputlabel = strcat('[',results_format,includebounds,']');
Results_securitiesmarket(fig_output,simtype,FIN_RESULTS,FIN_ASSETPRICES,T,T_sim,shocktime,FRFAtime,...
    results_format,includebounds,SM_outputlabel);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BANK FAILURES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results_failures(fig_output,simtype,n_banks,T_sim,FIN_FAILCOUNT,FIN_CUM_FAILS,...
    FIN_CAPITALSHORTFALL,FIN_NUM_ACTIVEBANKS,FIN_NUM_ACTIVEASSETS,...
    FIN_NUMEDGES,FIN_DENSITY,FIN_AVDEGREE,FIN_MU_A,FIN_MU_B,shocktime,FRFAtime);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
end
