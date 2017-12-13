function [banks,NMT_matrices,banksfail,n_banks,ibn_adjmat,opn_adjmat,...
    FailCount,FailedBankID,FailedBanks,ActiveBanks,num_ActiveBanks,num_FailedBanks,...
    ActiveAssets,num_ActiveAssets,num_InactiveAssets,CB_totalallotment] = ...
    Phase3(banks,DBV,DLV,NMT_matrices,banksfail,n_banks,m_assets,r_z,ibn_adjmat,opn_adjmat,AP,...
    FailCount,FailedBankID,FailedBanks,ActiveBanks,ActiveAssets,...
    CBintervention,refinancingfrequency,FRFA,CB_totalallotment,t,T,tol)

tau = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal policy intervention by central bank: Correct liquidity imbalances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = ActiveBanks
    banks(i).balancesheet.assets.cash(t,tau) = banks(i).balancesheet.assets.cash(t,tau-1);
    
   % if banks(i).balancesheet.liabilities.deposits(t,2) < 0
       % banks(i).balancesheet.liabilities.deposits(t,2) = banks(i).balancesheet.liabilities.deposits(t-1,2);
    %end
end

%------------------------------------------------
% Borrower refinancing
%------------------------------------------------

if strcmp(CBintervention,'on') && mod(t,refinancingfrequency) == 0
    
    loangap = zeros(1,numel(DBV));
    
    B_needliqcount = 0;
    
    for i = 1:numel(DBV)
        loangap(i) =  banks(DBV(i)).IBM.B_tot_requests(t) - banks(DBV(i)).IBM.B_tot_loans(t);
    
        if loangap(i) > tol
            B_needliqcount    = B_needliqcount+1;
            B_liqrec(B_needliqcount) = DBV(i);
        end
    end

    if B_needliqcount > 0
        B_CB_totalallotment = sum(loangap);

   
        for i = 1:numel(DBV)
            banks(DBV(i)).balancesheet.assets.cash(t,tau) = banks(DBV(i)).balancesheet.assets.cash(t,tau) +...
             loangap(i);
        end               
    else
        B_CB_totalallotment = 0;
    end

%------------------------------------------------
% Lender refinancing
%------------------------------------------------

   L_needliqcount = 0;
   repaymentgap = zeros(1,numel(DLV));

    for i = 1:numel(DLV)
        repaymentgap(i) =  banks(DLV(i)).IBM.L_tot_exp_repay(t) - banks(DLV(i)).IBM.L_tot_repaid_loans(t);
    
        if repaymentgap(i) > tol
            L_needliqcount    = L_needliqcount+1;
            L_liqrec(L_needliqcount) = DLV(i);
        end
    end

    if L_needliqcount > 0
        L_CB_totalallotment = sum(repaymentgap);
   
        for i = 1:numel(DLV)
            banks(DLV(i)).balancesheet.assets.cash(t,tau) = banks(DLV(i)).balancesheet.assets.cash(t,tau) +...
            repaymentgap(i);
        end               
    else
        L_CB_totalallotment = 0;
    end
elseif strcmp(CBintervention,'off') || mod(t,refinancingfrequency) ~= 0
    B_CB_totalallotment = 0;
    L_CB_totalallotment = 0;
end

CB_totalallotment(1) = B_CB_totalallotment;
CB_totalallotment(2) = L_CB_totalallotment;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st balance sheet update: Update total assets to account for returns on investment
% and central bank liquidity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = ActiveBanks

    % Assets after central bank intervention
        banks(i).balancesheet.assets.total(t,tau) = banks(i).balancesheet.assets.cash(t,tau)+...
        banks(i).balancesheet.assets.external_assets(t,tau);
    
    % Liabilities
    banks(i).balancesheet.liabilities.capital(t,tau) = banks(i).balancesheet.assets.total(t,tau)-...
        banks(i).balancesheet.liabilities.deposits(t,tau-2);
    
    banks(i).balancesheet.liabilities.total(t,tau) = banks(i).balancesheet.assets.total(t,tau);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine set of insolvent banks following CB intervention
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ea_holdings_mat = NMT_matrices(:,:,(4*t),1);
ea_port_mat     = NMT_matrices(:,:,(4*t),2);

for i = ActiveBanks
    
    banks(i).status(t) = '-';
    banks(i).failtime  = 0;
    
%     % Assets after end-of-period returns paid on investment
%     banks(i).balancesheet.assets.cash(t,tau) = banks(i).balancesheet.assets.cash(t,tau) +...
%             (r_z).*banks(i).balancesheet.assets.investment(t);
%         
%     banks(i).balancesheet.assets.total(t,tau) = banks(i).balancesheet.assets.cash(t,tau)+...
%         banks(i).balancesheet.assets.external_assets(t,tau);
    
    % Liabilities
    
   % banks(i).balancesheet.liabilities.capital(t,tau) = banks(i).balancesheet.assets.total(t,tau)-...
       % banks(i).balancesheet.liabilities.deposits(t,tau-2);
    
   % banks(i).balancesheet.liabilities.total(t,tau) = banks(i).balancesheet.assets.total(t,tau);
    
    if t<5
        insolvencyrequirement(i) =  false;
    else
        insolvencyrequirement(i) = logical(banks(i).balancesheet.liabilities.capital(t,tau)   < 0 &&...
                                    banks(i).balancesheet.liabilities.capital(t-1,tau) < 0 &&...
                                    banks(i).balancesheet.liabilities.capital(t-2,tau) < 0 &&...
                                    banks(i).balancesheet.liabilities.capital(t-3,tau) < 0 && ...
                                    banks(i).balancesheet.liabilities.capital(t-4,tau) < 0);
    end      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRFA policy active: CB intervenes by recapitalising insolvent banks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(FRFA,'on')
     for i = ActiveBanks
        
        if insolvencyrequirement(i)
            banks(i).balancesheet.assets.cash(t,tau) = banks(i).balancesheet.assets.cash(t,tau)+...
            abs(banks(i).balancesheet.liabilities.capital(t,tau))+...
            banks(i).balancesheet.liabilities.capital(1,1);
        
            FRFA_ind_allotment(i) = abs(banks(i).balancesheet.liabilities.capital(t,tau))+...
            banks(i).balancesheet.liabilities.capital(1,1);
        else
            FRFA_ind_allotment(i) = 0;
        end
        
        banks(i).status(t)   = 'A';
        banks(i).IBrolewhenfail = [];
        banks(i).failtime    = T;
        
     end
     
     CB_totalallotment(3) = sum(FRFA_ind_allotment);
     
     ibn_adjmat = ibn_adjmat;
     opn_adjmat = opn_adjmat;
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FRFA policy inactive: Insolvent banks are allowed to fail. External assets are
% divided amongst the bank's counterparties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif strcmp(FRFA,'off')
    
    CB_totalallotment(3) = 0;
    
    for i = ActiveBanks
        
        if insolvencyrequirement(i)
        
            % Banks for which insolvency condition is met: Report as failed
            FailCount    = FailCount + 1;
            FailedBankID(i) = 1;
        
            banks(i).status(t)      = 'F';
            banks(i).IBrolewhenfail = banks(i).IBM.status(t);
            banks(i).failtime       = t;
        
            fields = fieldnames(banks(i));
        
            for j = 1:length(fieldnames(banks))
                banksfail(i).(fields{j}) = banks(i).(fields{j});
            end
        
           % Divide external assets equally across failed bank' counterparties     
            for j = banks(i).counterpartyids
            
                banks(j).balancesheet.assets.external_assets(t,tau) = ...
                banks(j).balancesheet.assets.external_assets(t,tau)+...
                (1/banks(i).num_counterparties(t)).*(banks(i).balancesheet.assets.external_assets(t,tau));
            
               banks(j).balancesheet.assets.external_asset_holdings((4*t),:) = ones(1,banks(j).balancesheet.assets.num_external_assets(t)).*...
               (banks(j).balancesheet.assets.external_assets(t,tau)./banks(j).balancesheet.assets.num_external_assets(t));
           
                banks(j).balancesheet.assets.external_asset_port((4*t),:) = banks(j).balancesheet.assets.external_asset_holdings((4*t),:).*...
                AP(banks(j).balancesheet.assets.external_asset_ids);
            
                ea_holdings_mat(j,banks(j).balancesheet.assets.external_asset_ids) = banks(j).balancesheet.assets.external_asset_holdings(4*t,:);
                ea_port_mat(j,banks(j).balancesheet.assets.external_asset_ids)     = banks(j).balancesheet.assets.external_asset_port(4*t,:);
            end
            
            % Adjust external asset matrices of failed banks after redistribution across counterpartied.
            ea_holdings_mat(i,banks(i).balancesheet.assets.external_asset_ids) = banks(i).balancesheet.assets.num_external_assets(t);
            ea_port_mat(i,banks(i).balancesheet.assets.external_asset_ids)     = banks(i).balancesheet.assets.num_external_assets(t);
            
            banks(i).balancesheet.assets.external_assets(t,tau)                = 0;
                
            % Adjust adjacency matrices for interbank and overlapping portfolios to remove all connections of the failed bank
        
            ibn_adjmat(i,:) = zeros(1,n_banks);
            ibn_adjmat(:,i) = zeros(1,n_banks)';
        
            opn_adjmat(i,:) = zeros(1,m_assets);
        
        else 
            banks(i).status(t)      = 'A';
            banks(i).IBrolewhenfail = [];
            banks(i).failtime       = T;       
        end  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update set of active assets and banks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify assets that are no longer active (no holdings by any active banks)
    InactiveAssets = find(sum(opn_adjmat(:,ActiveAssets))==0);
    ActiveAssets   = setdiff(ActiveAssets,InactiveAssets);

    num_ActiveAssets   = numel(ActiveAssets);
    num_InactiveAssets = numel(InactiveAssets);

% Update set of active and failed banks based on current iteration of phase 3
    FailedBanks = unique([FailedBanks find(FailedBankID)]);
    ActiveBanks = setdiff(ActiveBanks,FailedBanks);

    num_ActiveBanks = numel(ActiveBanks);
    num_FailedBanks = numel(FailedBanks);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd balance sheet update:
%%% FRFA on: Liquidity injections for insolvent banks. Update cash
%%% FRFA off: Banks allowed to fail. External assets redistributed. Update external assets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = ActiveBanks
    
    banks(i).balancesheet.assets.cash(t,tau) = banks(i).balancesheet.assets.cash(t,tau) +...
    (r_z).*banks(i).balancesheet.assets.investment(t);

    % Assets after central bank intervention
    banks(i).balancesheet.assets.total(t,tau) = banks(i).balancesheet.assets.cash(t,tau)+...
    banks(i).balancesheet.assets.external_assets(t,tau);
    
    % Liabilities
    banks(i).balancesheet.liabilities.capital(t,tau) = banks(i).balancesheet.assets.total(t,tau)-...
        banks(i).balancesheet.liabilities.deposits(t,tau-2);
    
    banks(i).balancesheet.liabilities.total(t,tau) = banks(i).balancesheet.assets.total(t,tau);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update primary matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NMT_matrices(:,:,(4*t),1) = ea_holdings_mat(:,:);
NMT_matrices(:,:,(4*t),2) = ea_port_mat(:,:);

end