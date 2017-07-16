function [TOT_vec,TOT_sum,dTOT_vec,total_FS_vec,d_AP,d_FS,FailedBanks,ActiveBanks] =...
    Periodsummary(banks,ActiveBanks,FailedBanks,n_banks,IBratemat,TOT_vec,TOT_sum,dTOT_vec,...
    m_assets,assetprices,total_FS_vec,d_AP,d_FS,FailCount,FailedBankID,DBV,DLV,t,fileID_S)


%-------------------------------------------------------------------------
%% Initialisation
%-------------------------------------------------------------------------

total_assets_vec         = TOT_vec(:,1);
total_cash_vec           = TOT_vec(:,2);
total_ext_assets_vec     = TOT_vec(:,3);
total_des_investment_vec = TOT_vec(:,4);
total_investment_vec     = TOT_vec(:,5);

total_deposits_vec       = TOT_vec(:,6);
total_capital_vec        = TOT_vec(:,7);

total_requests_vec      = TOT_vec(:,8);
total_hoarding_vec      = TOT_vec(:,9);
total_loans_vec         = TOT_vec(:,10);
total_exp_repay_vec     = TOT_vec(:,11);
total_repayment_vec     = TOT_vec(:,12);

total_des_FS_vec        = TOT_vec(:,13);

total_assets         = TOT_sum(1,:);
total_cash           = TOT_sum(2,:);
total_ext_assets     = TOT_sum(3,:);

total_des_investment = TOT_sum(4,:);
total_investment     = TOT_sum(5,:);


total_deposits       = TOT_sum(6,:);
total_capital        = TOT_sum(7,:);

total_requests       = TOT_sum(8,:);
total_hoarding       = TOT_sum(9,:);
total_loans          = TOT_sum(10,:);
total_exp_repay      = TOT_sum(11,:);
total_repayment      = TOT_sum(12,:);

total_des_FS         = TOT_sum(13,:);
total_firesales      = TOT_sum(14,:);

d_total_assets       = dTOT_vec(1);
d_total_cash         = dTOT_vec(2);
d_total_investment   = dTOT_vec(3);
d_total_ext_assets   = dTOT_vec(4);
    
d_total_deposits     = dTOT_vec(5);
d_total_capital      = dTOT_vec(6);

%-------------------------------------------------------------------------
%% BALANCE SHEET OUTPUT
%-------------------------------------------------------------------------

fprintf(fileID_S,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\r\n');
fprintf(fileID_S,'--------------------------------------------------------------------- Period %d Summary  -----------------------------------------------------------------------\r\n',t);
fprintf(fileID_S,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\r\n');
fprintf(fileID_S,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');
fprintf(fileID_S,'Balance sheet\r\n');
fprintf(fileID_S,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');

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

total_assets(end)         = sum(total_assets_vec);
total_cash(end)           = sum(total_cash_vec);
total_ext_assets(end)     = sum(total_ext_assets_vec);
total_investment(end)     = nansum(total_investment_vec);
total_des_investment(end) = nansum(total_des_investment_vec);

total_deposits(end)       = sum(total_deposits_vec);
total_capital(end)        = sum(total_capital_vec);

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
    total_hoarding_vec(i)  = banks(i).IBM.hoarding(t);
    total_loans_vec(i)     = banks(i).IBM.L_tot_loans(t);
    total_exp_repay_vec(i) = banks(i).IBM.L_tot_exp_repay(t);
    total_repayment_vec(i) = banks(i).IBM.B_fin_tot_loanrepay(t);
       
end

total_requests(end)  = nansum(total_requests_vec);
total_hoarding(end)  = nansum(total_hoarding_vec);
total_loans(end)     = nansum(total_loans_vec);
total_exp_repay(end) = nansum(total_exp_repay_vec);
total_repayment(end) = nansum(total_repayment_vec);


fprintf(fileID_S,'Total interbank loan requests over %d borrowers = %.3f in period %d\r\n',numel(DBV),total_requests(end),t);
fprintf(fileID_S,'Total interbank loans (%.3f) + interest over %d lenders = %.3f\r\n',total_loans(end),numel(DLV),total_exp_repay(end));
fprintf(fileID_S,'Total loan repayment = %.3f\r\n',total_repayment(end));
    fprintf(fileID_S,'--> %.3f percent of total requests granted\r\n',(total_loans(end)/total_requests(end))*100);
    fprintf(fileID_S,'--> Lender hoarding constituted %.3f percent of total requests\r\n',(total_hoarding(end)/total_requests(end))*100);
    fprintf(fileID_S,'--> %.3f percent of total loans repaid\r\n',(total_repayment(end)/total_exp_repay(end))*100);
   
fprintf(fileID_S,'==================================================================================\r\n');
fprintf(fileID_S,'Rates\r\n');
fprintf(fileID_S,'==================================================================================\r\n');

    
IBratevec_temp = reshape(IBratemat,[1,n_banks^2]);
IBratevec = IBratevec_temp(IBratevec_temp~=0);

fprintf(fileID_S,'The average interbank rate was %.2f percent in period %d\r\n',mean(IBratevec),t);
fprintf(fileID_S,'[min max] = [%.2f %.2f]\r\n',min(IBratevec),max(IBratevec));

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
    total_FS_vec(i,:) = banks(i).firesales.ea_vec(t,:); 
    total_des_FS_vec(i)  = banks(i).firesales.tot_des_FS(t);    
end

asset_firesales      = sum(total_FS_vec);

total_des_FS(end)    = nansum(total_des_FS_vec);
total_firesales(end) = nansum(asset_firesales);

for i=1:m_assets   
    fprintf(fileID_S,'%.3f of asset %d sold at firesale in period %d\r\n',asset_firesales(i),i,t);    
end

if t>1    
    d_FS = ((total_firesales(end)-total_firesales(end-1))/total_firesales(end-1))*100;
else
    d_FS = 0;
end   

fprintf(fileID_S,'-----------------------------------------------------------------------------\r\n');     

fprintf(fileID_S,'Desired firesales =  %.3f\r\n',total_firesales(end));
fprintf(fileID_S,'Final firesales = %.3f\r\n',total_firesales(end));

if t>1
    fprintf(fileID_S,'--> Percent change = %.3f percent over last period\r\n',d_FS);
end

%-------------------------------------------------------------------------
%% BANK FAILURES
%-------------------------------------------------------------------------

failedbanks = sprintf(' %.3g ',find(FailedBankID));

fprintf(fileID_S,'%d banks failed in period %d\r\n', FailCount,t);
fprintf(fileID_S,'--> Failed banks IDs: [%s]\r\n',failedbanks);

FailedBanks = unique([FailedBanks find(FailedBankID)]);
ActiveBanks = setdiff(ActiveBanks,FailedBanks);

%-------------------------------------------------------------------------
%% Collecting output
%-------------------------------------------------------------------------

TOT_vec(:,1)  = total_assets_vec;
TOT_vec(:,2)  = total_cash_vec; 
TOT_vec(:,3)  = total_ext_assets_vec; 

TOT_vec(:,4)  = total_des_investment_vec;
TOT_vec(:,5)  = total_investment_vec;

TOT_vec(:,6)  = total_deposits_vec;
TOT_vec(:,7)  = total_capital_vec;

TOT_vec(:,8)  = total_requests_vec;
TOT_vec(:,9)  = total_hoarding_vec; 
TOT_vec(:,10) = total_loans_vec; 
TOT_vec(:,11) = total_exp_repay_vec;
TOT_vec(:,12) = total_repayment_vec;

TOT_vec(:,13) = total_des_FS_vec; 

TOT_sum(1,end)  = total_assets(end);
TOT_sum(2,end)  = total_cash(end); 
TOT_sum(3,end)  = total_ext_assets(end); 

TOT_sum(4,end)  = total_des_investment(end);
TOT_sum(5,end)  = total_investment(end);

TOT_sum(6,end)  = total_deposits(end);
TOT_sum(7,end)  = total_capital(end);

TOT_sum(8,end)  = total_requests(end);
TOT_sum(9,end)  = total_hoarding(end); 
TOT_sum(10,end) = total_loans(end); 
TOT_sum(11,end) = total_exp_repay(end);
TOT_sum(12,end) = total_repayment(end); 

TOT_sum(13,end) = total_des_FS(end);
TOT_sum(14,end) = total_firesales(end); 

dTOT_vec(1) = d_total_assets;
dTOT_vec(2) = d_total_cash;
dTOT_vec(3) = d_total_investment;
dTOT_vec(4) = d_total_ext_assets;
    
dTOT_vec(5) = d_total_deposits;
dTOT_vec(6) = d_total_capital;


end