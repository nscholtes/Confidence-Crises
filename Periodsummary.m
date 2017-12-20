function [Results_banks,Results_agg,Results_av,Results_min,Results_max,Results_dAgg,asset_FS_vec,d_AP,d_FS,num_edges,mu_A,mu_B] =...
    Periodsummary(banks,ActiveBanks,num_ActiveBanks,num_ActiveAssets,n_banks,IBratemat,...
    Results_banks,Results_agg,Results_av,Results_min,Results_max,Results_dAgg,...
    m_assets,assetprices,asset_FS_vec,d_AP,d_FS,ibn_adjmat,opn_adjmat,FailCount,FailedBankID,DBV,DLV,t,fileID_S)


%-------------------------------------------------------------------------
%% Initialisation
%-------------------------------------------------------------------------

% Individual bank information
total_assets_vec         = Results_banks(:,1);
total_cash_vec           = Results_banks(:,2);
total_ext_assets_vec     = Results_banks(:,3);
total_des_investment_vec = Results_banks(:,4);
total_investment_vec     = Results_banks(:,5);

total_deposits_vec       = Results_banks(:,6);
total_capital_vec        = Results_banks(:,7);

total_requests_vec      = Results_banks(:,8);
total_loans_BH_vec      = Results_banks(:,9);
total_hoarding_vec      = Results_banks(:,10);
total_loans_vec         = Results_banks(:,11);
total_exp_repay_vec     = Results_banks(:,12);
total_repayment_vec     = Results_banks(:,13);

total_des_FS_vec        = Results_banks(:,14);
total_act_FS_vec        = Results_banks(:,15);

% Aggregate across banks
total_assets         = Results_agg(1,:);
total_cash           = Results_agg(2,:);
total_ext_assets     = Results_agg(3,:);
total_des_investment = Results_agg(4,:);
total_investment     = Results_agg(5,:);

total_deposits       = Results_agg(6,:);
total_capital        = Results_agg(7,:);

total_requests      = Results_agg(8,:);
total_loans_BH      = Results_agg(9,:);
total_hoarding      = Results_agg(10,:);
total_loans         = Results_agg(11,:);
total_exp_repay     = Results_agg(12,:);
total_repayment     = Results_agg(13,:);

total_des_FS        = Results_agg(14,:);
total_act_FS        = Results_agg(15,:);

% Average over banks
av_assets         = Results_av(1,:);
av_cash           = Results_av(2,:);
av_ext_assets     = Results_av(3,:);
av_des_investment = Results_av(4,:);
av_investment     = Results_av(5,:);

av_deposits       = Results_av(6,:);
av_capital        = Results_av(7,:);

av_requests      = Results_av(8,:);
av_loans_BH      = Results_av(9,:);
av_hoarding      = Results_av(10,:);
av_loans         = Results_av(11,:);
av_exp_repay     = Results_av(12,:);
av_repayment     = Results_av(13,:);

av_des_FS        = Results_av(14,:);
av_act_FS        = Results_av(15,:);

av_IBrate        = Results_av(16,:);

% Min over banks
min_assets         = Results_min(1,:);
min_cash           = Results_min(2,:);
min_ext_assets     = Results_min(3,:);

min_des_investment = Results_min(4,:);
min_investment     = Results_min(5,:);

min_deposits       = Results_min(6,:);
min_capital        = Results_min(7,:);

min_requests       = Results_min(8,:);
min_loans_BH       = Results_min(9,:);
min_hoarding       = Results_min(10,:);
min_loans          = Results_min(11,:);
min_exp_repay      = Results_min(12,:);
min_repayment      = Results_min(13,:);

min_des_FS         = Results_min(14,:);
min_act_FS         = Results_min(15,:);

min_IBrate         = Results_min(16,:);


% Max over banks
max_assets         = Results_max(1,:);
max_cash           = Results_max(2,:);
max_ext_assets     = Results_max(3,:);

max_des_investment = Results_max(4,:);
max_investment     = Results_max(5,:);

max_deposits       = Results_max(6,:);
max_capital        = Results_max(7,:);

max_requests       = Results_max(8,:);
max_loans_BH       = Results_max(9,:);
max_hoarding       = Results_max(10,:);
max_loans          = Results_max(11,:);
max_exp_repay      = Results_max(12,:);
max_repayment      = Results_max(13,:);

max_des_FS         = Results_max(14,:);
max_act_FS         = Results_max(15,:);

max_IBrate         = Results_max(16,:);

% Change in aggregate variables
d_total_assets       = Results_dAgg(1);
d_total_cash         = Results_dAgg(2);
d_total_investment   = Results_dAgg(3);
d_total_ext_assets   = Results_dAgg(4);
    
d_total_deposits     = Results_dAgg(5);
d_total_capital      = Results_dAgg(6);

%-------------------------------------------------------------------------
%% BALANCE SHEET OUTPUT
%-------------------------------------------------------------------------

fprintf(fileID_S,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\r\n');
fprintf(fileID_S,'--------------------------------------------------------------------- Period %d Summary  -----------------------------------------------------------------------\r\n',t);
fprintf(fileID_S,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\r\n');
fprintf(fileID_S,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');
fprintf(fileID_S,'Balance sheet\r\n');
fprintf(fileID_S,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');

% Extract vectors from bank structures
for i = ActiveBanks
        
    % Assets  
    total_assets_vec(i)         = banks(i).balancesheet.assets.total(t,4);
    total_cash_vec(i)           = banks(i).balancesheet.assets.cash(t,4);
    total_ext_assets_vec(i)     = banks(i).balancesheet.assets.external_assets(t,4);
    total_des_investment_vec(i) = banks(i).balancesheet.assets.des_investment(t);
    total_investment_vec(i)     = banks(i).balancesheet.assets.investment(t);
    
    % Liabilities
    total_deposits_vec(i)       = banks(i).balancesheet.liabilities.deposits(t,2);
    total_capital_vec(i)        = banks(i).balancesheet.liabilities.capital(t,4);
               
end

% Summing over banks
total_assets(end)         = nansum(total_assets_vec);
total_cash(end)           = nansum(total_cash_vec);
total_ext_assets(end)     = nansum(total_ext_assets_vec);
total_investment(end)     = nansum(total_investment_vec);
total_des_investment(end) = nansum(total_des_investment_vec);

total_deposits(end)       = nansum(total_deposits_vec);
total_capital(end)        = nansum(total_capital_vec);

% Mean across banks
av_assets(end)         = nanmean(total_assets_vec); 
av_cash(end)           = nanmean(total_cash_vec);
av_ext_assets(end)     = nanmean(total_ext_assets_vec);
av_investment(end)     = nanmean(total_investment_vec);
av_des_investment(end) = nanmean(total_des_investment_vec);

av_deposits(end)       = nanmean(total_deposits_vec);
av_capital(end)        = nanmean(total_capital_vec);

% Min across banks
min_assets(end)         = min(total_assets_vec); 
min_cash(end)           = min(total_cash_vec);
min_ext_assets(end)     = min(total_ext_assets_vec);
min_investment(end)     = min(total_investment_vec);
min_des_investment(end) = min(total_des_investment_vec);

min_deposits(end)       = min(total_deposits_vec);
min_capital(end)        = min(total_capital_vec);

% Max across banks
max_assets(end)         = max(total_assets_vec); 
max_cash(end)           = max(total_cash_vec);
max_ext_assets(end)     = max(total_ext_assets_vec);
max_investment(end)     = max(total_investment_vec);
max_des_investment(end) = max(total_des_investment_vec);

max_deposits(end)       = max(total_deposits_vec);
max_capital(end)        = max(total_capital_vec);

if t>1
    
    d_total_assets     = ((total_assets(end) - total_assets(end-1))/(total_assets(end-1)))*100;
    d_total_cash       = ((total_cash(end) - total_cash(end-1))/(total_cash(end-1)))*100;
    d_total_ext_assets = ((total_ext_assets(end) - total_ext_assets(end-1))/(total_ext_assets(end-1)))*100;
    d_total_investment = ((total_investment(end) - total_investment(end-1))/(total_investment(end-1)))*100;
    
    d_total_deposits   = ((total_deposits(end) - total_deposits(end-1))/(total_deposits(end-1)))*100;
    d_total_capital    = ((total_capital(end) - total_capital(end-1))/(total_capital(end-1)))*100;
    
end

fprintf(fileID_S,'==================================================================================\r\n');
fprintf(fileID_S,'Assets\r\n');
fprintf(fileID_S,'==================================================================================\r\n');

fprintf(fileID_S,'Total assets = %.3f in period %d\r\n',total_assets(end),t);
    fprintf(fileID_S,'--> Represents a %.3f percent change over last period\r\n',d_total_assets);
    fprintf(fileID_S,'--> [min mean max] = [%.3f %.3f %.3f]\r\n', min(total_assets_vec),mean(total_assets_vec),max(total_assets_vec));
    
fprintf(fileID_S,'Total cash = %.3f in period %d\r\n',total_cash(end),t);
if t>1
    fprintf(fileID_S,'--> Represents a %.3f percent change over last period\r\n',d_total_cash);
end
    fprintf(fileID_S,'--> [min mean max] = [%.3f %.3f %.3f]\r\n', min(total_cash_vec),mean(total_cash_vec),max(total_cash_vec));
    
fprintf(fileID_S,'Total investments = %.3f in period %d\r\n',total_investment(end),t);
if t>1
    fprintf(fileID_S,'--> Represents a %.3f percent change over last period\r\n',d_total_investment);
end
    fprintf(fileID_S,'--> Total investment/Total desired investment = %.3f percent\r\n',(total_investment/total_des_investment)*100);
    
fprintf(fileID_S,'Total external assets = %.3f in period %d\n',total_ext_assets(end),t);
if t>1
    fprintf(fileID_S,'--> Represents a %.3f percent change over last period\r\n',d_total_ext_assets);
end
    fprintf(fileID_S,'--> [min mean max] = [%.3f %.3f %.3f]\r\n', min(total_ext_assets_vec),mean(total_ext_assets_vec),max(total_ext_assets_vec));
      
fprintf(fileID_S,'==================================================================================\r\n');
fprintf(fileID_S,'Liabilities\r\n');
fprintf(fileID_S,'==================================================================================\r\n');

fprintf(fileID_S,'Total deposits = %.3f in period %d\r\n',total_deposits(end),t);
if t>1
    fprintf(fileID_S,'--> Represents a %.3f percent change over last period\r\n',d_total_deposits);
end
    fprintf(fileID_S,'--> [min mean max] = [%.3f %.3f %.3f]\r\n', min(total_deposits_vec),mean(total_deposits_vec),max(total_deposits_vec));
    
fprintf(fileID_S,'Total capital = %.3f in period %d\r\n',total_capital(end),t);
if t>1
    fprintf(fileID_S,'--> Represents a %.3f percent change over last period\r\n',d_total_capital);
end
    fprintf(fileID_S,'--> [min mean max] = [%.3f %.3f %.3f]\r\n', min(total_capital_vec),mean(total_capital_vec),max(total_capital_vec));

fprintf(fileID_S,'-----------------------------------------------------------------------------\r\n');       

%-------------------------------------------------------------------------
%% INTERBANK MARKET OUTPUT
%-------------------------------------------------------------------------

fprintf(fileID_S,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');
fprintf(fileID_S,'Interbank market\r\n');
fprintf(fileID_S,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');

fprintf(fileID_S,'==================================================================================\r\n');
fprintf(fileID_S,'Volumes\r\n');
fprintf(fileID_S,'==================================================================================\r\n');


for i = ActiveBanks
            
    total_requests_vec(i)  = banks(i).IBM.L_tot_requests(t);
    total_loans_BH_vec(i)  = banks(i).IBM.L_tot_loans_BH(t);
    total_hoarding_vec(i)  = banks(i).IBM.hoarding(t);
    total_loans_vec(i)     = banks(i).IBM.L_tot_loans(t);
    total_exp_repay_vec(i) = banks(i).IBM.L_tot_exp_repay(t);
    total_repayment_vec(i) = banks(i).IBM.L_tot_repaid_loans(t);
    %total_repayment_vec(i) = banks(i).IBM.B_fin_tot_loanrepay(t);
       
end

IBrate_vec = nonzeros(IBratemat);

% Summing over banks
total_requests(end)  = nansum(total_requests_vec);
total_loans_BH(end)  = nansum(total_loans_BH_vec);
total_hoarding(end)  = nansum(total_hoarding_vec);
total_loans(end)     = nansum(total_loans_vec);
total_exp_repay(end) = nansum(total_exp_repay_vec);
total_repayment(end) = nansum(total_repayment_vec);

% Mean across banks
av_requests(end)  = nanmean(total_requests_vec);
av_loans_BH(end)  = nanmean(total_loans_BH_vec);
av_hoarding(end)  = nanmean(total_hoarding_vec);
av_loans(end)     = nanmean(total_loans_vec);
av_exp_repay(end) = nanmean(total_exp_repay_vec);
av_repayment(end) = nanmean(total_repayment_vec);

av_IBrate(end)  = mean(IBrate_vec);

% Min across banks
min_requests(end)  = min(total_requests_vec);
min_loans_BH(end)  = min(total_loans_BH_vec);
min_hoarding(end)  = min(total_hoarding_vec);
min_loans(end)     = min(total_loans_vec);
min_exp_repay(end) = min(total_exp_repay_vec);
min_repayment(end) = min(total_repayment_vec);

temp_min_IBrate = min(IBrate_vec);

if isempty(temp_min_IBrate)
    min_IBrate(end) = NaN;
else
    min_IBrate(end) = temp_min_IBrate;
end

% Max across banks
max_requests(end)  = max(total_requests_vec);
max_loans_BH(end)  = max(total_loans_BH_vec);
max_hoarding(end)  = max(total_hoarding_vec);
max_loans(end)     = max(total_loans_vec);
max_exp_repay(end) = max(total_exp_repay_vec);
max_repayment(end) = max(total_repayment_vec);

temp_max_IBrate = max(IBrate_vec);

if isempty(temp_max_IBrate)
    max_IBrate(end) = NaN;
else
    max_IBrate(end) = temp_max_IBrate;
end

fprintf(fileID_S,'Total interbank loan requests over %d borrowers = %.3f in period %d\r\n',numel(DBV),total_requests(end),t);
fprintf(fileID_S,'Total interbank loans (%.3f) + interest over %d lenders = %.3f\r\n',total_loans(end),numel(DLV),total_exp_repay(end));
fprintf(fileID_S,'Total loan repayment = %.3f\r\n',total_repayment(end));
    fprintf(fileID_S,'--> %.3f percent of total requests granted\r\n',(total_loans(end)/total_requests(end))*100);
    fprintf(fileID_S,'--> Lender hoarding constituted %.3f percent of total requests\r\n',(total_hoarding(end)/total_requests(end))*100);
    fprintf(fileID_S,'--> %.3f percent of total loans repaid\r\n',(total_repayment(end)/total_exp_repay(end))*100);
   
fprintf(fileID_S,'==================================================================================\r\n');
fprintf(fileID_S,'Rates\r\n');
fprintf(fileID_S,'==================================================================================\r\n');

%IBratevec_temp = reshape(IBratemat,[1,n_banks^2]);
%IBratevec = IBratevec_temp(IBratevec_temp~=0);

fprintf(fileID_S,'The average interbank rate was %.2f percent in period %d\r\n',av_IBrate,t);
fprintf(fileID_S,'[min max] = [%.2f %.2f]\r\n',min_IBrate,max_IBrate);

%-------------------------------------------------------------------------
%% SECURITIES MARKET
%-------------------------------------------------------------------------

fprintf(fileID_S,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');
fprintf(fileID_S,'Securities market\r\n');
fprintf(fileID_S,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');

fprintf(fileID_S,'==================================================================================\r\n');
fprintf(fileID_S,'Asset prices\r\n');
fprintf(fileID_S,'==================================================================================\r\n');

d_AP(1,:) = ((assetprices((4*t)-2,:)-assetprices((4*t)-3,:))./(assetprices((4*t)-3,:)))*100; % Due to APS
d_AP(2,:) = ((assetprices((4*t),:)-assetprices((4*t)-1,:))./(assetprices((4*t)-1,:)))*100; % Due to MIF

fprintf(fileID_S,'The average percent change in asset prices due to exogenous shocks = %.3f percent in period %d\r\n',mean(d_AP(1,:)),t);
fprintf(fileID_S,'The average percent change in asset prices due to 2nd round effects = %.3f percent in period %d\r\n',mean(d_AP(2,:)),t);

fprintf(fileID_S,'==================================================================================\r\n');
fprintf(fileID_S,'Firesales');
fprintf(fileID_S,'==================================================================================\r\n');

for i = ActiveBanks
    asset_FS_vec(i,:)    = banks(i).firesales.ea_vec(t,:);
    total_des_FS_vec(i)  = banks(i).firesales.tot_des_FS(t);    
    total_act_FS_vec(i)  = banks(i).firesales.final_firesales(t); 
end

asset_firesales   = sum(asset_FS_vec);

% Summing over banks
total_des_FS(end) = nansum(total_des_FS_vec);
total_act_FS(end) = nansum(total_act_FS_vec);

% Mean across banks
av_des_FS(end) = nanmean(total_des_FS_vec);
av_act_FS(end) = nanmean(total_act_FS_vec);

% Min across banks
min_des_FS(end) = min(total_des_FS_vec);
min_act_FS(end) = min(total_act_FS_vec);

% Max across banks
max_des_FS(end) = max(total_des_FS_vec);
max_act_FS(end) = max(total_act_FS_vec);

for k=1:m_assets   
    fprintf(fileID_S,'%.3f of asset %d sold at firesale in period %d\r\n',asset_firesales(k),k,t);    
end

if t>1    
    d_FS = ((total_act_FS(end)-total_act_FS(end-1))/total_act_FS(end-1))*100;
else
    d_FS = 0;
end   

fprintf(fileID_S,'-----------------------------------------------------------------------------\r\n');     

fprintf(fileID_S,'Desired firesales =  %.3f\r\n',total_act_FS(end));
fprintf(fileID_S,'Final firesales = %.3f\r\n',total_act_FS(end));

if t>1
    fprintf(fileID_S,'--> Percent change = %.3f percent over last period\r\n',d_FS);
end

%-------------------------------------------------------------------------
%% BANK FAILURES
%-------------------------------------------------------------------------

failedbanks = sprintf(' %.3g ',find(FailedBankID));

fprintf(fileID_S,'%d banks failed in period %d\r\n', FailCount,t);
fprintf(fileID_S,'--> Failed banks IDs: [%s]\r\n',failedbanks);

%-------------------------------------------------------------------------
%% NETWORK MEASURES
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interbank exposures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_edges = sum(sum(ibn_adjmat));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlapping portfolio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mu_A = Average degree of ASSETS
mu_A_vec = sum(opn_adjmat); 
mu_A     = sum(mu_A_vec)/num_ActiveAssets;
    
% mu_B = Average diversification of BANKS
mu_B_vec = sum(opn_adjmat,2);
mu_B     = sum(mu_B_vec)/num_ActiveBanks;
    
%-------------------------------------------------------------------------
%% Collecting output
%-------------------------------------------------------------------------

% TOT_vec: NxT matrices containing information for each bank over simulation
Results_banks(:,1)  = total_assets_vec;
Results_banks(:,2)  = total_cash_vec; 
Results_banks(:,3)  = total_ext_assets_vec; 

Results_banks(:,4)  = total_des_investment_vec;
Results_banks(:,5)  = total_investment_vec;

Results_banks(:,6)  = total_deposits_vec;
Results_banks(:,7)  = total_capital_vec;

Results_banks(:,8)  = total_requests_vec;
Results_banks(:,9)  = total_loans_BH_vec;
Results_banks(:,10)  = total_hoarding_vec; 
Results_banks(:,11) = total_loans_vec; 
Results_banks(:,12) = total_exp_repay_vec;
Results_banks(:,13) = total_repayment_vec;

Results_banks(:,14) = total_des_FS_vec; 
Results_banks(:,15) = total_act_FS_vec; 

% TOT_sum: 1xT matrices summing over all banks
Results_agg(1,end)  = total_assets(end);
Results_agg(2,end)  = total_cash(end); 
Results_agg(3,end)  = total_ext_assets(end); 

Results_agg(4,end)  = total_des_investment(end);
Results_agg(5,end)  = total_investment(end);

Results_agg(6,end)  = total_deposits(end);
Results_agg(7,end)  = total_capital(end);

Results_agg(8,end)  = total_requests(end);
Results_agg(9,end)  = total_loans_BH(end);
Results_agg(10,end)  = total_hoarding(end); 
Results_agg(11,end) = total_loans(end); 
Results_agg(12,end) = total_exp_repay(end);
Results_agg(13,end) = total_repayment(end); 

Results_agg(14,end) = total_des_FS(end);
Results_agg(15,end) = total_act_FS(end); 

% TOT_av: 1xT matrix containing average across banks in each time period
Results_av(1,end)  = av_assets(end);
Results_av(2,end)  = av_cash(end); 
Results_av(3,end)  = av_ext_assets(end); 

Results_av(4,end)  = av_des_investment(end);
Results_av(5,end)  = av_investment(end);

Results_av(6,end)  = av_deposits(end);
Results_av(7,end)  = av_capital(end);

Results_av(8,end)  = av_requests(end);
Results_av(9,end)  = av_loans_BH(end);
Results_av(10,end) = av_hoarding(end); 
Results_av(11,end) = av_loans(end); 
Results_av(12,end) = av_exp_repay(end);
Results_av(13,end) = av_repayment(end); 

Results_av(14,end) = av_des_FS(end);
Results_av(15,end) = av_act_FS(end); 

Results_av(16,end) = av_IBrate(end);

% TOT_min: 1xT matrix containing minimum across banks in each time period
Results_min(1,end)  = min_assets(end);
Results_min(2,end)  = min_cash(end); 
Results_min(3,end)  = min_ext_assets(end); 

Results_min(4,end)  = min_des_investment(end);
Results_min(5,end)  = min_investment(end);

Results_min(6,end)  = min_deposits(end);
Results_min(7,end)  = min_capital(end);

Results_min(8,end)  = min_requests(end);
Results_min(9,end)  = min_loans_BH(end);
Results_min(10,end) = min_hoarding(end); 
Results_min(11,end) = min_loans(end); 
Results_min(12,end) = min_exp_repay(end);
Results_min(13,end) = min_repayment(end); 

Results_min(14,end) = min_des_FS(end);
Results_min(15,end) = min_act_FS(end); 

Results_min(16,end) = min_IBrate(end);

% TOT_max: 1xT matrix containing minimum across banks in each time period
Results_max(1,end)  = max_assets(end);
Results_max(2,end)  = max_cash(end); 
Results_max(3,end)  = max_ext_assets(end); 

Results_max(4,end)  = max_des_investment(end);
Results_max(5,end)  = max_investment(end);

Results_max(6,end)  = max_deposits(end);
Results_max(7,end)  = max_capital(end);

Results_max(8,end)  = max_requests(end);
Results_max(9,end)  = max_loans_BH(end);
Results_max(10,end) = max_hoarding(end); 
Results_max(11,end) = max_loans(end); 
Results_max(12,end) = max_exp_repay(end);
Results_max(13,end) = max_repayment(end); 

Results_max(14,end) = max_des_FS(end);
Results_max(15,end) = max_act_FS(end); 

Results_max(16,end) = max_IBrate(end);

% Percent changes in variables between consecutive time periods
Results_dAgg(1) = d_total_assets;
Results_dAgg(2) = d_total_cash;
Results_dAgg(3) = d_total_investment;
Results_dAgg(4) = d_total_ext_assets;
    
Results_dAgg(5) = d_total_deposits;
Results_dAgg(6) = d_total_capital;

end