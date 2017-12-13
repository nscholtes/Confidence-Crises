function [ActiveBanks,ActiveAssets,banks,banksfail,num_ActiveBanks,num_FailedBanks,...
    num_ActiveAssets,num_InactiveAssets,FailCount_vec,FailedBankID_mat,FailedBanks,...
    numedges,mu_A,mu_B,assetprices,NT_matrices,NNT_matrices,NMT_matrices,TM_matrices,...
    IBN_adjmat,OPN_adjmat,Results_banks,Results_agg,Results_av,Results_min,Results_max,Results_dAgg,...
    total_firesales_vec,d_assetprices,d_firesales,CB_TOTallotment] =....
    siminitialize(T_sim,n_banks,m_assets,ibn_adjmat_init,opn_adjmat_init,IBN_adjmat,OPN_adjmat)

T = T_sim;

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

% Summarising and storing information after each round of the ABM
    Results_banks = zeros(n_banks,T,14); % Collecting balance sheet and interbank market information for each period
    Results_agg   = zeros(15,T);         % Aggregating across banks
    Results_av    = zeros(15,T);
    Results_min   = zeros(15,T);
    Results_max   = zeros(15,T);
    Results_dAgg  = zeros(6,T);

% Variables associated to bank firesales
    total_firesales_vec = zeros(n_banks,m_assets,T);
    d_assetprices       = zeros(2*T,m_assets);
    d_firesales         = zeros(1,T);
    
% Central bank allotment
    CB_TOTallotment = zeros(3,T);
    
end