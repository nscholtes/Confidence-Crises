function [banks,banksfail,n_banks,ibn_adjmat,opn_adjmat,FailCount,FailedBankID,ActiveAssets,num_IA,AP] = ...
    Phase3(banks,banksfail,n_banks,m_assets,r_z,ibn_adjmat,opn_adjmat,AP,...
    FailCount,FailedBankID,ActiveBanks,ActiveAssets,num_IA,t)

tau = 4;

for i = ActiveBanks
    
    banks(i).status(t) = '-';
    banks(i).failtime  = 0;
    
    % Assets after end-of-period returns paid on investment
    banks(i).balancesheet.assets.cash(t,tau) = banks(i).balancesheet.assets.cash(t,tau-1) +...
            (r_z).*banks(i).balancesheet.assets.investment(t);
        
    banks(i).balancesheet.assets.total(t,tau) = banks(i).balancesheet.assets.cash(t,tau)+...
        banks(i).balancesheet.assets.external_assets(t,tau);
    
    % Liabilities
    
    banks(i).balancesheet.liabilities.capital(t,tau) = banks(i).balancesheet.assets.total(t,tau)-...
        banks(i).balancesheet.liabilities.deposits(t,tau-2);
    
    banks(i).balancesheet.liabilities.total(t,tau) = banks(i).balancesheet.assets.total(t,tau);
        
end

% 1. Identify banks with negative end-of-period capital 

for i = ActiveBanks
        
    if banks(i).balancesheet.liabilities.capital(t,tau) < 0
        
        % Report bank as failed
        FailCount    = FailCount + 1;
        FailedBankID(i) = 1;
        
        banks(i).status(t)   = 'F';
        banks(i).IBrolewhenfail = banks(i).IBM.status(t);
        banks(i).failtime = t;
        
        fields = fieldnames(banks(i));
        
        for j = 1:length(fieldnames(banks))
            banksfail(i).(fields{j}) = banks(i).(fields{j});
        end
        
 % 2.  Divide external assets equally across failed bank' counterparties
        
        for j = banks(i).counterpartyids
            
            banks(j).balancesheet.assets.external_assets(t,tau) = ...
                banks(j).balancesheet.assets.external_assets(t,tau)+...
                (1/banks(i).num_counterparties(t)).*(banks(i).balancesheet.assets.external_assets(t,tau));
            
           banks(j).balancesheet.assets.external_asset_holdings((4*t),:) = ones(1,banks(j).balancesheet.assets.num_external_assets(t)).*...
               (banks(j).balancesheet.assets.external_assets(t,tau)./banks(j).balancesheet.assets.num_external_assets(t));
           
            banks(j).balancesheet.assets.external_asset_port((4*t),:) = banks(j).balancesheet.assets.external_asset_holdings((4*t),:).*...
                AP(banks(j).balancesheet.assets.external_asset_ids);
           
    
        end
        
        banks(i).balancesheet.assets.external_assets(t,tau) = 0;
                
% 3. Adjust adjacency matrices for interbank and overlapping portfolios to remove all connections of the failed bank
        
        ibn_adjmat(i,:) = zeros(1,n_banks);
        ibn_adjmat(:,i) = zeros(1,n_banks)';
        
        opn_adjmat(i,:) = zeros(1,m_assets);
        
    else 
        banks(i).status(t)   = 'A';
        banks(i).IBrolewhenfail = [];
        banks(i).failtime    = [];
        
        
        continue       
    end   
end

% Identify assets that are no longer active (no holdings by any active banks)

InactiveAssets = find(sum(opn_adjmat(:,ActiveAssets))==0);

if isempty(InactiveAssets)
    num_IA = 0;
else
    num_IA = numel(InactiveAssets);
end

% if t>1
%     IACount(2) = num_InactiveAssets - IACount(1);
% else
%     IACount(2) = num_InactiveAssets;
% end

%IACount = IACount(2);

% for i = ActiveBanks
%     for j = InactiveAssets
%         if ismember(j,banks(i).balancesheet.assets.external_asset_ids)
%             
%             Inactive_index = find(banks(i).balancesheet.assets.external_asset_ids==j)
%             banks(i).balancesheet.assets.external_asset_ids(Inactive_index) = [];
%             
%             banks(i).balancesheet.assets.external_asset_holdings((4*t),:) = ea_holdings_mat(i,banks(i).balancesheet.assets.external_asset_ids,2);
%             banks(i).balancesheet.assets.external_asset_port((4*t),:)     = ea_port_mat(i,banks(i).balancesheet.assets.external_asset_ids,2);
%             
%             AP(Inactive_index) = NaN;
%             
%         end
%     end
% end
            
ActiveAssets = setdiff(ActiveAssets,InactiveAssets);

%FailedBanks = unique([FailedBanks find(FailedBankID)]);
%ActiveBanks = setdiff(ActiveBanks,FailedBanks);

end