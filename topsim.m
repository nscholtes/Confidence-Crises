function [Networkdata, Graphdata,abmresults_allsims] = topsim()

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
%%% - Phase1.m: Encodes  phase 1 - Exogenous deposit shock and
%%%   move to interbank market
%%% - Phase2.m: Encodes phase 2 - Exogenous external asset shock, bank
%%%   firesales and loan repayment
%%% - Phase3.m: Encodes phase 3 - Bank insolvency computation and
%%%   central bank intervention

% Copyright (C) Nicolas K. Scholtes, 2017
% Distributed under GPL 3.0
%--------------------------------------------------------------------------
%% Script controls
%--------------------------------------------------------------------------
close all; clear; clc;

% Output figures to directory with LaTeX draft for automatic updating
fig_output = '/Users/nscholte/Desktop/Research/Ch.1 - Confidence crises/Drafts/Current version/Figures/';

% Create directory for simulation logs (with timestamp)

if ~exist('Logs','dir') && ~exist('Logs/Network Generation','dir')
    mkdir Logs
    mkdir Logs 'Network Generation' 
end

fileID_D   = fopen(strcat('Logs/detailed_log_',datestr(datetime('today')),'.txt'),'w');
fileID_S   = fopen(strcat('Logs/summary_log_',datestr(datetime('today')),'.txt'),'w');
fileID_IBN = fopen(strcat('Logs/Network Generation/IBN_',datestr(datetime('today')),'.txt'),'w');
%fileID_OPN = fopen(strcat('Logs/Network Generation/OPN_',datestr(datetime('today')),'.txt'),'w');

% Select whether to create a detailed simulation log  (Input: Y/N)
writeoption = 'N';

% Select whether to plot network in Gephi or MATLAB (Input: Gephi/MATLAB)
graphoption = 'MATLAB';

tol = 1.e-6; % Tolerance value for == conditional statements

timestamp = datestr(datetime('today'));

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
n_sims     = 5;     % Number of ABM simulations

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

num_shocked_ea_S = m_assets/5;
eashock_LB_S     = 0.9;
eashock_UB_S     = 1.1;
market_depth_S   = zeros(1,m_assets); % Market depth of assets used in firesale computation
    
num_shocked_ea_C = m_assets;
eashock_LB_C     = 0.75;
eashock_UB_C     = 1;
market_depth_C   = ones(1,m_assets);
    
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

simtype        = 'crisis';
MRR            = 0.02;  % Minimum reserve requirement
CBintervention = 'on';
FRFA           = 'off'; % Fixed rate full allotment: Banks have unlimited access to CB liquidity


if strcmp(FRFA,'on')
    FRFA_label = '_FRFAon';
elseif strcmp(FRFA,'off')
    FRFA_label = '_FRFAoff';
end

% Combining the simulation labels for later output

allsimlabels = strcat('[',simtype,FRFA_label,']');

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
    median_density = mean(density_vec(median_index));

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

IBN_adjmat = zeros(n_banks,n_banks,T);
OPN_adjmat  = zeros(n_banks,m_assets,T);

%--------------------------------------------------------------------------
% Initialisation
%--------------------------------------------------------------------------

% Summarising and storing information after each round of the ABM
Results_banks = zeros(n_banks,T,14); % Collecting balance sheet and interbank market information for each period
Results_agg   = zeros(15,T);         % Aggregating across banks
Results_av    = zeros(15,T);
Results_min   = zeros(15,T);
Results_max   = zeros(15,T);
Results_dAgg  = zeros(6,T);

% Storing summarised matrix for each run of the ABM
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

sim_CB_TOTallotment = zeros(n_sims,T,2);

close all 
clc
   
%-----------------------------------------------------------------------------------------------------------------------------------------------------
%% Agent-Based Model
%-----------------------------------------------------------------------------------------------------------------------------------------------------
tic
for k = 1:n_sims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialising/resetting output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define label to store cleaned ABM results for each run
    abmresults_label{k} = strcat('sim',num2str(k));

%  Variables associated to the set of banks and external assets
    ActiveBanks        = 1:n_banks;
    ActiveAssets       = 1:m_assets;
    banks              = struct;
    banksfail          = struct;
    num_ActiveBanks    = zeros(1,T);
    num_FailedBanks    = zeros(1,T);
    num_ActiveAssets   = zeros(1,T);
    num_InactiveAssets = zeros(1,T);
    FailCount_vec      = zeros(1,T);
    FailedBankID_mat   = zeros(T,n_banks);
    FailedBanks        = [];
    numedges           = zeros(1,T);
    mu_A               = zeros(1,T);
    mu_B               = zeros(1,T);

% Assetprice vector: 4 steps per time period: 1) initialisation, 2) asset price shock, 3) firesales, 4) market impact function
    assetprices = ones(4*T,m_assets);

% Matrices storing relevant information across time periods 
    NT_matrices     = zeros(n_banks,T,3);             % Shock distribution information
    NNT_matrices    = zeros(n_banks,n_banks,T,8);     % Matrix storing interbank variable dynamics
    NMT_matrices    = zeros(n_banks,m_assets,4*T,2);  % Matrix storing external asset dynamics
    TM_matrices     = zeros(T,m_assets,2);            % Matrix storing variables for market impact function following firesales 

% Map initial simulated adjacency matrix to all time periods (to be updated during the ABM)
    for t = 1:T
        IBN_adjmat(:,:,t) = ibn_adjmat_init;
        OPN_adjmat(:,:,t) = opn_adjmat_init; 
    end

% Variables associated to bank firesales
    total_firesales_vec = zeros(n_banks,m_assets,T);
    d_assetprices       = zeros(2*T,m_assets);
    d_firesales         = zeros(1,T);
    
% Central bank allotment
    CB_TOTallotment = zeros(2,T);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixing shock structure
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
        elseif strcmp(simtype,'baseline')
            Pars_pshock_current = Pars_pshock(1,:);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        disp('--------------- SIMULATION PARAMETERS ---------------')
        fprintf('The current simulation type is: %s\n',simtype)
        if strcmp(simtype,'crisis')
            fprintf('--> Crisis starts at period %d and lasts %d periods\n',shocktime(1),shocktime(end)-shocktime(1))
        end
        fprintf('Central bank intervention is %s\n',CBintervention)
        fprintf('--> Intervention type is: %s\n',FRFA)
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
            CBintervention,FRFA,CB_TOTallotment(:,t),t,T,tol); 
        
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
    end
    
    % Balance sheet and interbank market dynamics
    sim_Results_agg(:,:,k) = Results_agg;
    sim_Results_av(:,:,k)  = Results_av;
    sim_Results_min(:,:,k) = Results_min;
    sim_Results_max(:,:,k) = Results_max;
    sim_assetprices(:,:,k) = assetprices;
    
    % Interbank rate dynamics
    
    
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
    
    abmresults_allsims.(abmresults_label{k}).vars = abmresults(banks,n_banks,MRR);
    abmresults_allsims.(abmresults_label{k}).ActiveBanks = ActiveBanks;
    abmresults_allsims.(abmresults_label{k}).FailedBanks = FailedBanks;
    
    abmresults_allsims.(abmresults_label{k}).ActiveAssets = ActiveAssets;

clearvars banks DBV DLV B_noIB 
clc
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COLLECTING SIMULUATION RESULTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fin_Results,fin_assetprices,fin_FailCount,fin_cum_Fails,fin_capitalshortfall,fin_num_ActiveBanks,fin_num_ActiveAssets,...
    fin_numnodes,fin_numedges,fin_density,fin_avdegree,fin_mu_A,fin_mu_B,fin_CB_TOTallotment] = ...
    simcollect(abmresults_allsims,abmresults_label,n_sims,sim_Results_agg,sim_Results_av,sim_Results_min,sim_Results_max,...
    sim_assetprices,sim_failcount,sim_num_ActiveBanks,sim_num_ActiveAssets,sim_numnodes,sim_numedges,...
    sim_density,sim_avdegree,sim_mu_A,sim_mu_B,T_sim,sim_CB_TOTallotment);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLEAN UP  AND SAVE WORKSPACE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
clearvars a a_max a_min alpha av_div beta d_assetprices d_firesales DBV DLV eashock_LB_C eashock_LB_S...
    eashock_UB_C eashock_UB_S fileID_D fileID_IBN fileID_S gamma dTOT_matrices market_depth_C...
    market_depth_S num_shocked_ea_C num_shocked_ea_S Pars_balancesheet Pars_dshock Pars_interestrates...
    Pars_netgen Pars_opnet Pars_policy Pars_pshock Pars_pshock_current psi r_b r_d r_e r_z SPF theta tol

save(strcat('output_',timestamp,allsimlabels,'.mat'));                            % Save workspace as .mat file

%--------------------------------------------------------------------------
%% Model output
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NETWORK VISUALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fix parameters used for network visualization
numslices    = 4;        % Number of time periods to use for network visualization
timeslices   = linspace(T/numslices,T,numslices); 
NNTvisualize = [1 5 8];  % Select which weighted matrices to visualize: (1 = ex-post adjacency matrix, 5 = Loan matrix, 8 = Interbank rate matrix)
NTvisualize  = 1;        % Total assets of active banks in each period
NMTvisualize = 2;
markernorm   = 10;       % Normalization factor such that marker size of largest node = $markernorm$
edgewnorm    = 10;       % Normalization factor such that width of largest edge = %edgewnorm$
graphlayout  = 'circle'; % Choose from:  'auto' (default) | 'circle' | 'force' | 'layered' | 'subspace' | 'force3' | 'subspace3'

networkparameters = [markernorm edgewnorm graphlayout];

[Networkdata,Graphdata] = viewnetworks(fig_output,networkparameters,Results_banks(:,[1 timeslices],NTvisualize),[1 timeslices],...
    cat(3,ibn_adjmat_init,IBN_adjmat(:,:,timeslices)),NNT_matrices(:,:,[1 timeslices],NNTvisualize),cat(3,opn_adjmat_init,OPN_adjmat(:,:,timeslices)));
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
Results_balancesheet(fig_output,fin_Results,...
    T,T_sim,normalizegraph,results_format,includebounds,BS_outputlabel,allsimlabels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERBANK MARKET DYNAMICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IBM_outputlabel = strcat('[',results_format,includebounds,']');
Results_IBM(fig_output,fin_Results,fin_CB_TOTallotment,T_sim,results_format,includebounds,IBM_outputlabel,allsimlabels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECURITIES MARKET DYNAMICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results_securitiesmarket(fig_output,fin_Results,fin_assetprices,T,allsimlabels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BANK FAILURES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results_failures(fig_output,banksfail,fin_FailCount,fin_cum_Fails,fin_capitalshortfall,fin_num_ActiveBanks,fin_num_ActiveAssets,...
    fin_numedges,fin_density,fin_avdegree,fin_mu_A,fin_mu_B,allsimlabels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
