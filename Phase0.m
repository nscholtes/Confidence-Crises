function [banks,NMT_matrices,assetprices,tau] =...
            Phase0(banks,ActiveBanks,a,BS_pars,T,t,ibn_adjmat,NMT_matrices,opn_adjmat,assetprices)

% -------------------------------------------------------------------------
%% Initialisation
% -------------------------------------------------------------------------

alpha = BS_pars(1);
beta  = BS_pars(2);

if t == 1
    ea_holdings_mat = NMT_matrices(:,:,t,1);   
    ea_port_mat     = NMT_matrices(:,:,t,2);   
else
    ea_holdings_mat = NMT_matrices(:,:,(4*t)-4:(4*t)-3,1);
    ea_port_mat     = NMT_matrices(:,:,(4*t)-4:(4*t)-3,2);
end

tau = 1;

% -------------------------------------------------------------------------
%% Inter-period parameter update
% -------------------------------------------------------------------------

% First period inputs are exogenous, determined by the network simulation stage.
% For t=2,...T, parameters updated within each period are carried over to the next.
 
if t == 1
        
    for i = ActiveBanks  
        
    % Specifying bank identifiers
    
        banks(i).id = i;
        
    % Balance sheet: Asset-side
        
        % Initialising
        banks(i).balancesheet.assets.total           = zeros(T,4);
        banks(i).balancesheet.assets.cash            = zeros(T,4);
        
        banks(i).balancesheet.assets.external_assets = zeros(T,4);

        % Populating the first entry from node size distribution
        banks(i).balancesheet.assets.total(t,tau)           = a(i);
        banks(i).balancesheet.assets.cash(t,tau)            = (1-alpha)*a(i);
        banks(i).balancesheet.assets.external_assets(t,tau) = alpha*a(i);
  
    % Balance sheet: Liabilities-side
        
        banks(i).balancesheet.liabilities.total    = zeros(T,4);
        banks(i).balancesheet.liabilities.deposits = zeros(T,2);
        banks(i).balancesheet.liabilities.capital  = zeros(T,4);
        
        banks(i).balancesheet.liabilities.total(t,tau)    = a(i);
        banks(i).balancesheet.liabilities.deposits(t,tau) = beta*a(i);
        banks(i).balancesheet.liabilities.capital(t,tau)  = (1-beta)*a(i);
        
    % Counterparty information
    
        banks(i).neighbour_vec(1,:,t)  = ibn_adjmat(i,:);
        banks(i).counterpartyids       = find(banks(i).neighbour_vec(1,:,1));
        banks(i).num_counterparties(t) = numel(banks(i).counterpartyids);
        
    % External asset portfolios (overlapping)
    
        banks(i).balancesheet.assets.external_asset_ids   = find(opn_adjmat(i,:)==1);  
        banks(i).balancesheet.assets.num_external_assets  = ...
            ones(1,T)*numel(banks(i).balancesheet.assets.external_asset_ids);

        banks(i).balancesheet.assets.external_asset_port(t,:) = ones(1,sum(opn_adjmat(i,:)))...
            .*(banks(i).balancesheet.assets.external_assets(t,tau))./sum(opn_adjmat(i,:));
        
        banks(i).balancesheet.assets.external_asset_holdings(t,:) = banks(i).balancesheet.assets.external_asset_port(t,:);
        
        for k = t+1:4*T
            banks(i).balancesheet.assets.external_asset_port(k,:) = zeros(1,numel(banks(i).balancesheet.assets.external_asset_port(t,:)));
            banks(i).balancesheet.assets.external_asset_holdings(k,:)  = zeros(1,numel(banks(i).balancesheet.assets.external_asset_holdings(t,:)));
        end
        
        
        ea_holdings_mat(i,banks(i).balancesheet.assets.external_asset_ids) = banks(i).balancesheet.assets.external_asset_port(t,:);
        ea_port_mat(i,banks(i).balancesheet.assets.external_asset_ids)     = ea_holdings_mat(i, banks(i).balancesheet.assets.external_asset_ids,t);
        
        banks(i).depositshock = zeros(1,T);
                        
    end
            
    NMT_matrices(:,:,1,1) = ea_holdings_mat;
    NMT_matrices(:,:,1,2) = ea_port_mat;
    
elseif t > 1
    
    assetprices((4*t)-3,:) = assetprices((4*t)-4,:);
    
    ea_holdings_mat(:,:,2) = ea_holdings_mat(:,:,1);
    ea_port_mat(:,:,2)     = ea_port_mat(:,:,1); 
         
    for i = ActiveBanks
        
% Counterparty information
    
        banks(i).neighbour_vec(1,:,t)  = ibn_adjmat(i,:);
        banks(i).counterpartyids       = find(banks(i).neighbour_vec(1,:,t));
        banks(i).num_counterparties(t) = numel(banks(i).counterpartyids);
        
% Asset-side
        
        banks(i).balancesheet.assets.cash(t,tau)            = banks(i).balancesheet.assets.cash(t-1,4);
        banks(i).balancesheet.assets.external_assets(t,tau) = banks(i).balancesheet.assets.external_assets(t-1,4);
        banks(i).balancesheet.assets.total(t,tau)           = banks(i).balancesheet.assets.total(t-1,4);
          
% External asset portfolio

        banks(i).balancesheet.assets.external_asset_ids      = find(opn_adjmat(i,:)==1); 
        banks(i).balancesheet.assets.num_external_assets(t)  = numel(banks(i).balancesheet.assets.external_asset_ids);

        banks(i).balancesheet.assets.external_asset_holdings((4*t)-3,:) = banks(i).balancesheet.assets.external_asset_holdings((4*t)-4,:);
        banks(i).balancesheet.assets.external_asset_port((4*t)-3,:)     = banks(i).balancesheet.assets.external_asset_port((4*t)-4,:);
                
% Liabilities-side 
    
        banks(i).balancesheet.liabilities.deposits(t,tau) = banks(i).balancesheet.liabilities.deposits(t-1,2);
        banks(i).balancesheet.liabilities.total(t,tau)    = banks(i).balancesheet.liabilities.total(t-1,4);
        
        banks(i).balancesheet.liabilities.capital(t,tau)  = banks(i).balancesheet.liabilities.capital(t-1,4);
  
    end
    
    NMT_matrices(:,:,(4*t)-3,1) = ea_holdings_mat(:,:,2);
    NMT_matrices(:,:,(4*t)-3,2) = ea_port_mat(:,:,2);
    
end

end
