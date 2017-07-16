function [banks,assetprices,opn_adjmat,NMT_matrices,NNT_matrices,TM_matrices]...
            =Phase2(banks,ActiveBanks,ActiveAssets,opn_adjmat,Pars_pshock_current,Pars_policy,...
            NMT_matrices,NNT_matrices,TM_matrices,assetprices,DBV,DLV,t,fileID_D,writeoption)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase 2: External asset shock drives down the market prices of banks'securities portfolios and affects their ability to sell
% off assets to meet interbank obligations from Phase 1 and comply with the calibrated regulatory constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
%% INITIALISATION
%-------------------------------------------------------------------------

tau = 3;

% Import state variables from topsim used in current phase

num_shocked_ea = Pars_pshock_current(1);
eashock_LB     = Pars_pshock_current(2);
eashock_UB     = Pars_pshock_current(3);
market_depth   = Pars_pshock_current(4);

MRR = Pars_policy;

% Import interbank matrices from the current period

borr_rep_mat        = NNT_matrices(:,:,t,4);
adj_lender_lend_mat = NNT_matrices(:,:,t,6);

% Import external asset matrices (initialised in phase 0 along with empty rows to be filled in current phase)

ea_holdings_mat = NMT_matrices(:,:,(4*t)-3:4*t,1); 
ea_port_mat     = NMT_matrices(:,:,(4*t)-3:4*t,2);

tot_eaFS_vec = TM_matrices(t,:,1); % Read in zero vectors
tot_eaH_vec  = TM_matrices(t,:,2);

for i = ActiveBanks
    ea_holdings_mat(i,:,2) = ea_holdings_mat(i,:,1);
    ea_holdings_mat(i,:,3) = ea_holdings_mat(i,:,2);
    ea_holdings_mat(i,:,4) = ea_holdings_mat(i,:,3);
    
    ea_port_mat(i,:,2)  = ea_port_mat(i,:,1);
    ea_port_mat(i,:,3)  = ea_port_mat(i,:,2);
    ea_port_mat(i,:,4)  = ea_port_mat(i,:,3);
    
    banks(i).firesales.ea_vec(t,ActiveAssets) = zeros(1,numel(ActiveAssets));
    banks(i).balancesheet.assets.eaH_vec(t,ActiveAssets) = zeros(1,numel(ActiveAssets));
end

%-------------------------------------------------------------------------
%% EXTERNAL ASSET SHOCK - FIRST ROUND EFFECTS
%------------------------------------------------------------------------- 

% Select subset of assets to shock - random integer draw

num_shocked_ea = numel(ActiveAssets);

for i = 1:num_shocked_ea
    shocked_ea_vec(i)= randi(numel(ActiveAssets));   
end

non_shocked_ea_vec = setdiff(linspace(1,numel(ActiveAssets),numel(ActiveAssets)),shocked_ea_vec);

% Apply the random shock to the selected subset of assets 

assetprices((4*t)-2,shocked_ea_vec)     = (eashock_LB+(eashock_UB -eashock_LB).*rand(1,num_shocked_ea)).*assetprices((4*t)-3,shocked_ea_vec);
assetprices((4*t)-2,non_shocked_ea_vec) = assetprices((4*t)-3,non_shocked_ea_vec);

for i = ActiveBanks  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% 1st Update: Asset price shock affects value of each bank' portfolio (but no change in holdings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    banks(i).shocked_external_assets = intersect(shocked_ea_vec,banks(i).balancesheet.assets.external_asset_ids);
    
    ea_port_mat(i,banks(i).shocked_external_assets,2) = ea_holdings_mat(i,banks(i).shocked_external_assets,2)...
        .*assetprices((4*t)-2,banks(i).shocked_external_assets);    
   
    % Updating external asset entries in bank structures
    banks(i).balancesheet.assets.external_asset_holdings((4*t)-2,:) = ea_holdings_mat(i,banks(i).balancesheet.assets.external_asset_ids,2);
    banks(i).balancesheet.assets.external_asset_port((4*t)-2,:)     = ea_port_mat(i,banks(i).balancesheet.assets.external_asset_ids,2);
    
    banks(i).balancesheet.assets.external_assets(t,tau-1) = sum(banks(i).balancesheet.assets.external_asset_port((4*t)-2,:));
 
 end

%-------------------------------------------------------------------------
%% LOAN REPAYMENT AND FIRESALES
%------------------------------------------------------------------------- 

myprint(writeoption,fileID_D,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');
myprint(writeoption,fileID_D,'Phase 2: Loan repayment and firesales\r\n');
myprint(writeoption,fileID_D,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% Borrowers determine how much liquidity they have available for repaying interbank loans
myprint(writeoption,fileID_D,'==================================================================================\r\n');
myprint(writeoption,fileID_D,'Step 1: Determine whether firesales are required\r\n');
myprint(writeoption,fileID_D,'==================================================================================\r\n');

assetprices((4*t)-1,:) = assetprices((4*t)-2,:); % No change in asset prices during firesales

for i = 1:numel(DBV)
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contingency for borrowers who did not have any lending counterparties in
% the current period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    banks(DBV(i)).firesales.final_firesales(t) = 0;
    
    banks(DBV(i)).IBM.L_bil_repaid_loans    = zeros(1,banks(DBV(i)).num_final_cps(t));
    banks(DBV(i)).IBM.L_tot_repaid_loans(t) = NaN;
    
    if isempty(banks(DBV(i)).final_cps) % Case for borrowers who did not have lending counterparties in the current period
        
        banks(DBV(i)).IBM.B_req_bil_loanrepay      = 0;
        banks(DBV(i)).IBM.B_fin_bil_loanrepay      = 0;
        
        banks(DBV(i)).IBM.B_req_tot_loanrepay(t)   = 0; 
        banks(DBV(i)).IBM.B_fin_tot_loanrepay(t)   = 0;
        banks(DBV(i)).IBM.B_prov_tot_loanrepay(t)  = 0;
      
        banks(DBV(i)).IBM.canrepay_NoFS(t)         = NaN;
        banks(DBV(i)).firesales.fullrepaywithFS(t) = NaN;
        banks(DBV(i)).firesales.act_FS_vec         = zeros(1,banks(DBV(i)).balancesheet.assets.num_external_assets(t));
        banks(DBV(i)).firesales.final_firesales(t) = 0;
        
        banks(DBV(i)).firesales.tot_des_FS(t)      = 0;
        
        banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,:) = ...
            banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-2,:);
                    
        banks(DBV(i)).balancesheet.assets.external_asset_port((4*t)-1,:) =...
            banks(DBV(i)).balancesheet.assets.external_asset_port((4*t)-2,:);
             
        banks(DBV(i)).balancesheet.assets.external_assets(t,tau) = sum(banks(DBV(i)).balancesheet.assets.external_asset_port((4*t)-1,:));
        
        banks(DBV(i)).balancesheet.assets.cash(t,tau) = banks(DBV(i)).balancesheet.assets.cash(t,tau-1);
       
        myprint(writeoption,fileID_D,'Borrower %d has no lenders in period %d: no firesales\r\n',DBV(i),t);
        myprint(writeoption,fileID_D,'----------------------------------------------------------------------------------\r\n');
            
    else
    
    % Borrowers compute bilateral and total loan obligations (obtained from current period loan+interest matrix)
    
        temp_repaystore = adj_lender_lend_mat(:,DBV(i))';
        temp_repaystore(temp_repaystore==0)=[];
    
        banks(DBV(i)).IBM.B_req_bil_loanrepay    = temp_repaystore;
        banks(DBV(i)).IBM.B_req_tot_loanrepay(t) = sum(banks(DBV(i)).IBM.B_req_bil_loanrepay);
    
        clearvars temp_repaystore
    
        myprint(writeoption,fileID_D,'Borrower %d has total loans of %.3f to repay across %d lenders in period %d\r\n',...
            DBV(i),banks(DBV(i)).IBM.B_req_tot_loanrepay(t),banks(DBV(i)).num_final_cps(t),t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scenario in which banks have enough liquidity to repay loan in full: No firesale of external assets needed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        if (banks(DBV(i)).balancesheet.assets.cash(t,tau-1))*(1-MRR) >= banks(DBV(i)).IBM.B_req_tot_loanrepay(t)
           
            banks(DBV(i)).IBM.canrepay_NoFS(t) = 1;
            banks(DBV(i)).firesales.fullrepaywithFS(t)  = NaN;
        
            banks(DBV(i)).IBM.B_fin_bil_loanrepay     = banks(DBV(i)).IBM.B_req_bil_loanrepay;
            banks(DBV(i)).IBM.B_fin_tot_loanrepay(t)  = banks(DBV(i)).IBM.B_req_tot_loanrepay(t);
            
            myprint(writeoption,fileID_D,'--> Borrower %d has sufficient reserves (%.3f) to repay loans in period %d\r\n',...
                DBV(i),(banks(DBV(i)).balancesheet.assets.cash(t,tau-1))*(1-MRR),t);
            
            % 1st end of phase update (BORROWERS): Borrowers pay off loans without resorting to firesales 

            banks(DBV(i)).balancesheet.assets.cash(t,tau) = banks(DBV(i)).balancesheet.assets.cash(t,tau-1)...
                - banks(DBV(i)).IBM.B_fin_tot_loanrepay(t);
            
            myprint(writeoption,fileID_D,'--> Updated reserves: %.3f\r\n',banks(DBV(i)).balancesheet.assets.cash(t,tau));
            myprint(writeoption,fileID_D,'----------------------------------------------------------------------------------\r\n');
            
            borr_rep_mat(DBV(i),banks(DBV(i)).final_cps) = banks(DBV(i)).IBM.B_fin_bil_loanrepay;
                    
            banks(DBV(i)).firesales.act_FS_vec = zeros(1,banks(DBV(i)).balancesheet.assets.num_external_assets(t));
        
            banks(DBV(i)).firesales.final_firesales(t) = 0;
            banks(DBV(i)).firesales.tot_des_FS(t)      = 0;
            
            % No change to external asset portfolio holdings and value (only market price changes)
            
            banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,:) = ...
                banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-2,:);
                    
             banks(DBV(i)).balancesheet.assets.external_asset_port((4*t)-1,:) =...
                 banks(DBV(i)).balancesheet.assets.external_asset_port((4*t)-2,:);
             
             banks(DBV(i)).balancesheet.assets.external_assets(t,tau) = sum(banks(DBV(i)).balancesheet.assets.external_asset_port((4*t)-1,:));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Inadequate liquidity: Borrower engages in firesale of external assets to repay loans
% Sell off assets at prevailing market price (post-shock)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Outer loop: Determine whether firesales are required

        elseif (banks(DBV(i)).balancesheet.assets.cash(t,tau-1))*(1-MRR)  <  banks(DBV(i)).IBM.B_req_tot_loanrepay(t)
        
            banks(DBV(i)).IBM.canrepay_NoFS(t) = 0;
        
            banks(DBV(i)).firesales.tot_des_FS(t) =  banks(DBV(i)).IBM.B_req_tot_loanrepay(t)-(banks(DBV(i)).balancesheet.assets.cash(t,tau-1))*(1-MRR);
        
            banks(DBV(i)).firesales.des_FS_vec    =  ones(1,banks(DBV(i)).balancesheet.assets.num_external_assets(t))...
                .*(1./banks(DBV(i)).balancesheet.assets.num_external_assets(t)).*banks(DBV(i)).firesales.tot_des_FS(t);
        
            banks(DBV(i)).firesales.act_FS_vec = zeros(1,banks(DBV(i)).balancesheet.assets.num_external_assets(t));
        
            myprint(writeoption,fileID_D,'Borrower %d has insufficient reserves to repay interbank loans (%.3f<%.3f). Has to firesale %.3f of each of its %d external assets\r\n',...
                DBV(i),...
                (banks(DBV(i)).balancesheet.assets.cash(t,tau-1))*(1-MRR),...
                banks(DBV(i)).IBM.B_req_tot_loanrepay(t),...
                (1./banks(DBV(i)).balancesheet.assets.num_external_assets(t)).*banks(DBV(i)).firesales.tot_des_FS(t),...
                banks(DBV(i)).balancesheet.assets.num_external_assets(t));
        
%-------------------------------------------------FIRESALES-----------------------------------------------------------

% STEP 1: Engage in firesales of individual external assets

% Inner loop: As banks sell off assets, check that holdings are sufficient. 
% Loop through each asset and sell off until bilateral repayment requirement is met
            
            for j = 1:banks(DBV(i)).balancesheet.assets.num_external_assets(t)
             
% When required firesale of an asset is > holding of that asset, sell off entire position of that asset:
        
                if  banks(DBV(i)).firesales.des_FS_vec(j)/assetprices((4*t)-1,banks(DBV(i)).balancesheet.assets.external_asset_ids(j))...
                        > banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-2,j)
                    
                    %disp('Case')
                
                    banks(DBV(i)).firesales.act_FS_vec(j) = banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,j);
                
                    banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,j) = 0;
                    
                    banks(DBV(i)).balancesheet.assets.external_asset_port((4*t)-1,j) = ...
                        banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,j).*...
                        assetprices((4*t)-1,banks(DBV(i)).balancesheet.assets.external_asset_ids(j));
                    
                    opn_adjmat(DBV(i),banks(DBV(i)).balancesheet.assets.external_asset_ids(j)) = 0;
                    
                    myprint(writeoption,fileID_D,'--> Sell off entire position of asset %d to meet firesale requirement\r\n',...
                        banks(DBV(i)).balancesheet.assets.external_asset_ids(j));
                
% Current asset holdings sufficient for firesales:
                
                elseif banks(DBV(i)).firesales.des_FS_vec(j)/assetprices((4*t)-1,banks(DBV(i)).balancesheet.assets.external_asset_ids(j))...
                        <= banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-2,j)
               
                    banks(DBV(i)).firesales.act_FS_vec(j) = banks(DBV(i)).firesales.des_FS_vec(j);
                    
                    banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,j) = banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-2,j)-...
                        banks(DBV(i)).firesales.act_FS_vec(j);  
                    
                    banks(DBV(i)).balancesheet.assets.external_asset_port((4*t)-1,j) =...
                        banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,j).*...
                        assetprices((4*t)-1,banks(DBV(i)).balancesheet.assets.external_asset_ids(j));
                    
                    myprint(writeoption,fileID_D,'--> Sufficient holdings (%.3f) of asset %d to meet requirement (%.3f). Updated holdings given by: %.3f\r\n',...
                        banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-2,j),...
                        banks(DBV(i)).balancesheet.assets.external_asset_ids(j),...
                        (1/banks(DBV(i)).balancesheet.assets.num_external_assets(t)).*banks(DBV(i)).firesales.tot_des_FS(t),...
                        banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,j));
                else
                    disp('Error in firesale mechanism!')
                end            
            end
            
            banks(DBV(i)).balancesheet.assets.external_assets(t,tau) = sum(banks(DBV(i)).balancesheet.assets.external_asset_port((4*t)-1,:));
                    
%---------------------------------------------------------------------------------------------------------------------
% STEP 2: Determine whether firesales allowed borrowers to satisfy interbank obligations
 
        banks(DBV(i)).firesales.final_firesales(t) = sum(banks(DBV(i)).firesales.act_FS_vec);
        banks(DBV(i)).IBM.B_prov_tot_loanrepay(t) = banks(DBV(i)).firesales.final_firesales(t) + (banks(DBV(i)).balancesheet.assets.cash(t,tau-1))*(1-MRR);
        
% Insufficient liquidity even with firesales: Repay interbank loans on a pro-rata basis
 
            if banks(DBV(i)).IBM.B_prov_tot_loanrepay(t) < banks(DBV(i)).IBM.B_req_tot_loanrepay(t) 
            
                banks(DBV(i)).firesales.fullrepaywithFS(t) = 0;
                banks(DBV(i)).IBM.B_fin_bil_loanrepay    = banks(DBV(i)).IBM.B_prov_tot_loanrepay(t)./banks(DBV(i)).IBM.B_req_tot_loanrepay(t)...
                .*banks(DBV(i)).IBM.B_req_bil_loanrepay;
                banks(DBV(i)).IBM.B_fin_tot_loanrepay(t) = sum(banks(DBV(i)).IBM.B_fin_bil_loanrepay);
             
                banks(DBV(i)).firesales.fullrepaywithFS(t) = NaN;
                
% With firesales, borrowers are able to fulfill their interbank obligations
                
            elseif banks(DBV(i)).IBM.B_prov_tot_loanrepay(t) >= banks(DBV(i)).IBM.B_req_tot_loanrepay(t) 
                
                banks(DBV(i)).IBM.B_fin_bil_loanrepay    = banks(DBV(i)).IBM.B_req_bil_loanrepay;
                banks(DBV(i)).IBM.B_fin_tot_loanrepay(t) = banks(DBV(i)).IBM.B_req_tot_loanrepay(t);
                
                banks(DBV(i)).firesales.fullrepaywithFS(t)  = 1;
                
            else
                disp('Error in borrower repayment with firesales!')
            end
              
        myprint(writeoption,fileID_D,'Borrower %d firesales %.3f worth of external assets in period %d to pay back %.3f worth of interbank loans. Repayment/loan ratio is: %.3f percent\r\n',...
        DBV(i),banks(DBV(i)).firesales.final_firesales(t),t,banks(DBV(i)).IBM.B_req_tot_loanrepay(t),...
        (banks(DBV(i)).IBM.B_fin_tot_loanrepay(t)/banks(DBV(i)).IBM.B_req_tot_loanrepay(t))*100);
    
       % 1st end of phase update (BORROWERS): Remaining reserves after loan repayment and firesales 
        banks(DBV(i)).balancesheet.assets.cash(t,tau) = banks(DBV(i)).balancesheet.assets.cash(t,tau-1)...
            -  banks(DBV(i)).IBM.B_fin_tot_loanrepay(t) + banks(DBV(i)).firesales.final_firesales(t);
        
        myprint(writeoption,fileID_D,'--> Updated reserves: %.3f\r\n',banks(DBV(i)).balancesheet.assets.cash(t,tau));        
        myprint(writeoption,fileID_D,'-----------------------------------------------------------------------------\r\n');  
        end
    end
       
    borr_rep_mat(DBV(i),banks(DBV(i)).final_cps) = banks(DBV(i)).IBM.B_fin_bil_loanrepay;
    
    banks(DBV(i)).IBM.L_bil_repaid_loans(t) = 0;
    banks(DBV(i)).IBM.L_tot_repaid_loans(t) = 0;
    
    banks(DBV(i)).firesales.external_asset_firesales(t,:) = banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-2,:)-...
        banks(DBV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,:);
    
   banks(DBV(i)).firesales.ea_vec(t,banks(DBV(i)).balancesheet.assets.external_asset_ids) = banks(DBV(i)).firesales.external_asset_firesales(t,:);
    
end

%-------------------------------------------------------------------------
%% UPDATING FIRESALE ENTRIES FOR LENDERS (TRIVIAL)
%------------------------------------------------------------------------- 

for i =1:numel(DLV)
    
banks(DLV(i)).balancesheet.assets.external_asset_holdings((4*t)-1,:) = banks(DLV(i)).balancesheet.assets.external_asset_holdings((4*t)-2,:);
    banks(DLV(i)).balancesheet.assets.external_asset_port((4*t)-1,:) = banks(DLV(i)).balancesheet.assets.external_asset_port((4*t)-2,:); 
    
    banks(DLV(i)).balancesheet.assets.external_assets(t,tau) = sum(banks(DLV(i)).balancesheet.assets.external_asset_port((4*t)-1,:));
    
    banks(DLV(i)).firesales.external_asset_firesales(t,:) = zeros(1,banks(DLV(i)).balancesheet.assets.num_external_assets(t));
        
end
%-------------------------------------------------------------------------
%% MARKET IMPACT FUNCTION - SECOND ROUND EFFECTS
%------------------------------------------------------------------------- 

% 1. Compute total firesales and total holdings (across all banks) for each individual asset in current period

for i = ActiveBanks   
    temp_tot_eaFS_vec(i,:) = banks(i).firesales.ea_vec(t,:);
    
    banks(i).balancesheet.assets.eaH_vec(t,banks(i).balancesheet.assets.external_asset_ids) =...
        banks(i).balancesheet.assets.external_asset_holdings((4*t)-1,:);
    
    temp_tot_eaH_vec(i,:) = banks(i).balancesheet.assets.eaH_vec(t,:);
end

tot_eaH_vec_1  = sum(temp_tot_eaH_vec);
tot_eaFS_vec_1 = sum(temp_tot_eaFS_vec);

clearvars temp_tot_eaFS_vec temp_tot_eaH_vec

keepassets = find(tot_eaH_vec_1);

tot_eaH_vec  = tot_eaH_vec_1(keepassets);
tot_eaFS_vec = tot_eaFS_vec_1(keepassets);

% 2. MIF: Exponential function of firesales

market_depth = ones(1,numel(keepassets));

assetprices((4*t),:) = assetprices((4*t)-1,:);

assetprices(4*t,keepassets) = assetprices((4*t)-1,keepassets).*exp(-(market_depth.*tot_eaFS_vec)./tot_eaH_vec);

for i = ActiveBanks
    
    banks(i).balancesheet.assets.external_asset_holdings((4*t),:) = banks(i).balancesheet.assets.external_asset_holdings((4*t)-1,:);
    banks(i).balancesheet.assets.external_asset_port((4*t),:)     = banks(i).balancesheet.assets.external_asset_holdings((4*t),:).*...
        assetprices(4*t,banks(i).balancesheet.assets.external_asset_ids);
  
    banks(i).balancesheet.assets.external_assets(t,tau+1) = sum(banks(i).balancesheet.assets.external_asset_port((4*t),:));
    
end

%-------------------------------------------------------------------------
%% LENDERS REACT TO BORROWER LOAN REPAYMENT
%------------------------------------------------------------------------- 

myprint(writeoption,fileID_D,'==================================================================================\r\n');
myprint(writeoption,fileID_D,'Step 3: Lenders react to borrower loan repayment\r\n');
myprint(writeoption,fileID_D,'==================================================================================\r\n');

% Updating external assets entry for LENDERS. Do not firesale assets (no change in holdings) but price shock changes balance sheet entry

for i = 1:numel(DLV)
    
    banks(DLV(i)).IBM.canrepay_NoFS(t)         = NaN;
    banks(DLV(i)).firesales.fullrepaywithFS(t) = NaN;  
    banks(DLV(i)).firesales.final_firesales(t) = NaN;
    banks(DLV(i)).firesales.tot_des_FS(t)      = NaN;
    
    banks(DLV(i)).IBM.B_req_tot_loanrepay(t)  = NaN;
    banks(DLV(i)).IBM.B_req_bil_loanrepay(t)  = NaN;

    banks(DLV(i)).IBM.B_prov_tot_loanrepay(t) = NaN;
    banks(DLV(i)).IBM.B_fin_tot_loanrepay(t)  = NaN;
    banks(DLV(i)).IBM.B_fin_bil_loanrepay(t)  = NaN;
   
    temp_repstore = borr_rep_mat(:,DLV(i))';
    temp_repstore(temp_repstore==0) = [];
    
    banks(DLV(i)).IBM.L_bil_repaid_loans  = temp_repstore;

    clearvars temp_repstore
  
    banks(DLV(i)).IBM.L_tot_repaid_loans(t)   = sum(banks(DLV(i)).IBM.L_bil_repaid_loans);
    
    % 1st end of phase update (LENDERS): Borrower loan repayment added to reserves
    
    banks(DLV(i)).balancesheet.assets.cash(t,tau) = banks(DLV(i)).balancesheet.assets.cash(t,tau-1) + ...
        banks(DLV(i)).IBM.L_tot_repaid_loans(t);
    
    if banks(DLV(i)).numborrowing_cps(t) == 0
        
        myprint(writeoption,fileID_D,'Lender %d had no borrowing counterparties in period %d\r\n',DLV(i),t);       
        myprint(writeoption,fileID_D,'-----------------------------------------------------------------------------\r\n');     
        
    else
    
        myprint(writeoption,fileID_D,'Lender %d lent out %.3f (= %.3f with interest) to %d counterparties and received %.3f in period %d\r\n',...
            DLV(i),banks(DLV(i)).IBM.L_tot_loans(t),banks(DLV(i)).IBM.L_tot_exp_repay(t),...
            banks(DLV(i)).numborrowing_cps(t),banks(DLV(i)).IBM.L_tot_repaid_loans(t),t);
    
        myprint(writeoption,fileID_D,'Lender %d total repayment to loan ratio = %.3f percent in period %d\r\n',...
            DLV(i),(banks(DLV(i)).IBM.L_tot_repaid_loans(t)./banks(DLV(i)).IBM.L_tot_exp_repay(t))*100,t);
    
    
        for j = 1: banks(DLV(i)).numborrowing_cps(t)
        
            myprint(writeoption,fileID_D,'--> Lender %d and Borrower %d loan/repayment ratio = %.3f percent\r\n',...
                DLV(i),banks(DLV(i)).borrowing_cps(j),(banks(DLV(i)).IBM.L_bil_repaid_loans(j)./banks(DLV(i)).IBM.L_bil_exp_repay(j))*100);        
        end
        
       myprint(writeoption,fileID_D,'-----------------------------------------------------------------------------\r\n');       
    end
end
%-------------------------------------------------------------------------
%% BALANCE SHEET AND HOLDING MATRICES UPDATING
%------------------------------------------------------------------------- 

for i = ActiveBanks

% External asset portfolio
    % 1st update: External asset shock
    ea_holdings_mat(i,banks(i).balancesheet.assets.external_asset_ids,2) = ...
        banks(i).balancesheet.assets.external_asset_holdings((4*t)-2,:);
    
    ea_port_mat(i,banks(i).balancesheet.assets.external_asset_ids,2)  = ...
        banks(i).balancesheet.assets.external_asset_port((4*t)-2,:);
    
    %2nd update: Asset firesales to meet interbank obligations
     ea_holdings_mat(i,banks(i).balancesheet.assets.external_asset_ids,3) = ...
        banks(i).balancesheet.assets.external_asset_holdings((4*t)-1,:);
    
    ea_port_mat(i,banks(i).balancesheet.assets.external_asset_ids,3)  = ...
        banks(i).balancesheet.assets.external_asset_port((4*t)-1,:);
    
    %3rd update: 2nd round price spirals via market impact function
     ea_holdings_mat(i,banks(i).balancesheet.assets.external_asset_ids,4) = ...
        banks(i).balancesheet.assets.external_asset_holdings((4*t),:);
    
    ea_port_mat(i,banks(i).balancesheet.assets.external_asset_ids,4)  = ...
        banks(i).balancesheet.assets.external_asset_port((4*t),:);

% Assets
    
    banks(i).balancesheet.assets.total(t,tau) = ...
          banks(i).balancesheet.assets.cash(t,tau)...
        + banks(i).balancesheet.assets.investment(t)...
        + banks(i).balancesheet.assets.external_assets(t,tau+1);
    
% Liabilities    
    
    banks(i).balancesheet.liabilities.capital(t,tau) =  banks(i).balancesheet.assets.total(t,tau)...
        -  banks(i).balancesheet.liabilities.deposits(t,tau-1);
    
    banks(i).balancesheet.liabilities.total(t,tau) = banks(i).balancesheet.assets.total(t,tau);
  
end

%-------------------------------------------------------------------------
%% COLLECTING OUTPUT
%-------------------------------------------------------------------------

NNT_matrices(:,:,t,4) = borr_rep_mat';
NNT_matrices(:,:,t,6) = adj_lender_lend_mat;

NMT_matrices(:,:,(4*t)-3:4*t,1) = ea_holdings_mat(:,:,1:4);
NMT_matrices(:,:,(4*t)-3:4*t,2) = ea_port_mat(:,:,1:4);

TM_matrices(t,keepassets,1) = tot_eaFS_vec; 
TM_matrices(t,keepassets,2) = tot_eaH_vec;

end
