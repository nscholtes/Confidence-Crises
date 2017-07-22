function [Networkdata, Graphdata] = topsim()

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

%--------------------------------------------------------------------------
% CALIBRATION OF SIMULATION PARAMETERS
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bank size/fitness distribution parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_banks     = 50;        % Number of banks
ActiveBanks = 1:n_banks;
a_min       = 5;         % Minimum bank size
a_max       = 100;       % Maximum bank size
gamma       = 2;         % Power law exponent for bank size distribution

Pars_netgen = [a_min,a_max,gamma,n_banks];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlapping portfolio parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_assets = 50; % Number of external assets
av_div   = 10;   % Average number of assets held per bank (diversification)

ActiveAssets = 1:m_assets;

Pars_opnet = [m_assets, av_div];

%--------------------------------------------------------------------------
% Agent-Based Model Parameters
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model timing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 1000;        % Number of iteration steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial balance sheet weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.9;   % external asset/asset ratio
beta  = 0.92;  % deposit/liability ratio

Pars_balancesheet= [alpha beta];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Liquidity and portfolio shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

information = 'perfect';

% Liquidity (deposit shock)

theta = 0.5;   % Multiplier for mean-reverting component
SPF   = 0.025; % Shock proportionality factor

% Portfolio shock
shocktime   = T/2; 

num_shocked_ea_S = m_assets/5;
eashock_LB_S     = 0.9;
eashock_UB_S     = 1.1;
market_depth_S   = 0.2*ones(1,m_assets); % Market depth of assets used in firesale computation
    
num_shocked_ea_C = m_assets;
eashock_LB_C     = 0.75;
eashock_UB_C     = 1;
market_depth_C    = 0.5*ones(1,m_assets);
    
%Collecting parameters
Pars_dshock  = [theta SPF];
Pars_pshock = [num_shocked_ea_S, eashock_LB_S, eashock_UB_S market_depth_S;
    num_shocked_ea_C, eashock_LB_C, eashock_UB_C market_depth_C];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interest rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_z  = 0.05; % Illiquid investment rate
r_d  = 0.02; % Deposit rate
r_e  = 0.03; % External asset return
r_b  = 0.04; % Initial interbank rat

Pars_interestrates = [r_z r_d r_e r_b];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Policy parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Macroprudential policy parameters + central bank instruments

%CR = 0.08;
%LR = 0.03;

psi = 0.02;  % Minimum reserve requirement

Pars_policy = psi;

%--------------------------------------------------------------------------
% Generating the network
%--------------------------------------------------------------------------

[a,ibn_adjmat_init] = netgen(n_banks,gamma,a_min,a_max,fileID_IBN,fig_output);
[opn_adjmat_init]   = opnet(n_banks, m_assets, av_div,T,fig_output);

IBN_adjmat = zeros(n_banks,n_banks,T);
OPN_adjmat  = zeros(n_banks,m_assets,T);

for t = 1:T
    IBN_adjmat(:,:,t) = ibn_adjmat_init;
    OPN_adjmat(:,:,t) = opn_adjmat_init;  
end

%--------------------------------------------------------------------------
% Initialisation
%--------------------------------------------------------------------------

banks           = struct;
banksfail       = struct;
assetprices     = ones(4*T,m_assets);
NT_matrices     = zeros(n_banks,T,3);             % Shock distribution information
NNT_matrices    = zeros(n_banks,n_banks,T,8);     % Matrix storing interbank variable dynamics
NMT_matrices    = zeros (n_banks,m_assets,4*T,2); % Matrix storing external asset dynamics
TM_matrices     = zeros(T,m_assets,2);            % Matrix storing variables for market impact function following firesales 
TOT_NT_matrices = zeros(n_banks,T,13);            % Collecting balance sheet and interbank market information for each period
TOT_T_matrices  = zeros(14,T);
dTOT_matrices   = zeros(6,T);

total_firesales_vec = zeros(n_banks,m_assets,T);
d_assetprices       = zeros(2*T,m_assets);
d_firesales         = zeros(1,T);

FailCount_vec    = zeros(1,T);
FailedBankID_mat = zeros(T,n_banks);
FailedBanks      = [];

num_InactiveAssets = zeros(1,T);

close all 
clc
   
%--------------------------------------------------------------------------
%% Agent-Based Model
%--------------------------------------------------------------------------
tic
for t=1:T
    if t==1
        x = 1;
    else  
        x = t-1;
    end        
    if t < shocktime
        Pars_pshock_current = Pars_pshock(1,:);
    elseif t >= shocktime
        Pars_pshock_current = Pars_pshock(2,:);
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODEL DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(1,'Current period: %d\n',t);
% Phase 0: Initialisation, Population of the network and carrying over of variables from previous iterations
   [banks,NMT_matrices,assetprices] =...
        Phase0(banks,ActiveBanks,a,Pars_balancesheet,T,t,IBN_adjmat(:,:,x),NMT_matrices,OPN_adjmat(:,:,x),assetprices);
% Phase 1: Deposit shock, investment decision, interbank market
   [banks,DBV,DLV,NT_matrices,NNT_matrices] = ...
        Phase1(banks,ActiveBanks,IBN_adjmat(:,:,x),n_banks,beta,Pars_interestrates,Pars_dshock,Pars_policy,information,NT_matrices,NNT_matrices,t,fileID_D,writeoption);
% Phase 2: Asset price shock, loan repayment, firesales and second round effects
   [banks,assetprices,OPN_adjmat(:,:,t),NMT_matrices,NNT_matrices,TM_matrices] = ...
        Phase2(banks,ActiveBanks,ActiveAssets,OPN_adjmat(:,:,x),Pars_pshock_current,Pars_policy,...
        NMT_matrices,NNT_matrices,TM_matrices,assetprices,DBV,DLV,t,fileID_D,writeoption);
% Phase 3: Removal of insolvent banks and (possible) central bank intervention
   [banks,banksfail,n_banks,IBN_adjmat(:,:,t),OPN_adjmat(:,:,t),FailCount_vec(t),FailedBankID_mat(t,:),ActiveAssets,num_InactiveAssets(t),assetprices(t,:)] =  ...
       Phase3(banks,banksfail,n_banks,m_assets,Pars_interestrates(1),IBN_adjmat(:,:,x),OPN_adjmat(:,:,x),assetprices(t,:),...
       FailCount_vec(t),FailedBankID_mat(t,:),ActiveBanks,ActiveAssets,num_InactiveAssets(t),t,T);      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CURRENT ITERATION SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [TOT_NT_matrices(:,t,:),TOT_T_matrices(:,1:t),dTOT_matrices(:,t),...
        total_firesales_vec(:,:,t),d_assetprices((2*t)-1:2*t,:),d_firesales(t),FailedBanks,ActiveBanks] =...
        Periodsummary(banks,ActiveBanks,FailedBanks,n_banks,NNT_matrices(:,:,t,8),TOT_NT_matrices(:,t,:),TOT_T_matrices(:,1:t),dTOT_matrices(:,t),...
        m_assets,assetprices,total_firesales_vec(:,:,t),d_assetprices((2*t)-1:2*t,:),d_firesales(t),...
        FailCount_vec(t),FailedBankID_mat(t,:),DBV,DLV,t,fileID_S);    
    if numel(ActiveBanks) == 0
        fprintf(1,'All banks failed by period %d\n',t)
        T = t;
        break
    end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTRA-PERIOD HOUSEKEEPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    [banks] = housekeeping(banks,n_banks);
    clc;
end
toc

abmresults = abmresults(banks,n_banks);
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

[Networkdata,Graphdata] = viewnetworks(fig_output,networkparameters,TOT_NT_matrices(:,[1 timeslices],NTvisualize),[1 timeslices],...
    cat(3,ibn_adjmat_init,IBN_adjmat(:,:,timeslices)),NNT_matrices(:,:,[1 timeslices],NNTvisualize),cat(3,opn_adjmat_init,OPN_adjmat(:,:,timeslices)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATION OF TABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%allparameters = [a_min,a_max,gamma,av_div,BS_parsvec,calibraterate_vec,n_banks,m_assets];
%maketables(allparameters);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION RESULTS TO COMMAND WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%allparameters = [a_min,a_max,gamma,av_div,BS_parsvec,calibraterate_vec,n_banks,m_assets];
%allvariables  = [ActiveBanks,ActiveAssets,TOT_T_matrices,IBN_adjmat(:,:,[1,T])];

%simoutput(allparameters,allvariables);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BALANCE SHEET DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalizegraph = 'Y'; % Y/N to normalize balance sheet variables by total assets
Results_balancesheet(fig_output,TOT_T_matrices,TOT_NT_matrices,normalizegraph)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTERBANK MARKET DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results_IBM(fig_output,TOT_T_matrices,TOT_NT_matrices,NNT_matrices(:,:,:,8),T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECURITIES MARKET DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Results_securitiesmarket(fig_output,TOT_T_matrices,TM_matrices,total_firesales_vec,assetprices,T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BANK FAILURE AND CENTRAL BANK INTERVENTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[banksfail] = Results_failures(fig_output,n_banks,IBN_adjmat,OPN_adjmat,banksfail,FailedBanks,FailCount_vec,num_InactiveAssets,t,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars a a_max a_min alpha av_div beta d_assetprices d_firesales...
DBV DLV dTOT_matrices eashock_LB_C eashock_LB_S eashock_UB_C eashock_UB_S gamma m_assets market_depth_C market_depth_S...
MRR num_shocked_ea_C num_shocked_ea_S r_b r_d r_e r_z shocktime SPF statevar_vec theta

end
