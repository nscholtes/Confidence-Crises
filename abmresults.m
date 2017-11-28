function[abmresults] = abmresults(banks,n_banks,MRR)

LA_variables  = {'Number of borrowers';'Total requests';'Available reserves';'Provisional loans';...
    'Hoarding';'Total loans';'Expected repayment';'Final repayment';'Reserves after repayment'};
BA_variables  = {'Number of lenders';'Total requests';'Total loans';'Required repayment';'Available reserves';...
    'Desired firesales';'Final firesales';'Final repayment';'Reserves after repayment'};
BS_variables  = {'Cash';'External assets';'Deposits';'Capital'};
Inv_variables = {'Desired investment';'Final investment'};

for i = 1:n_banks
% Balance sheet
    cash_array = [banks(i).balancesheet.assets.cash(1,1) banks(i).balancesheet.assets.cash(1:banks(i).failtime,4)'];
    ea_array   = [banks(i).balancesheet.assets.external_assets(1,1) banks(i).balancesheet.assets.external_assets(1:banks(i).failtime,4)'];
    
    dep_array  = [banks(i).balancesheet.liabilities.deposits(1,1) banks(i).balancesheet.liabilities.deposits(1:banks(i).failtime,2)'];
    cap_array  = [banks(i).balancesheet.liabilities.capital(1,1) banks(i).balancesheet.liabilities.capital(1:banks(i).failtime,4)'];

    BS_array = [cash_array; ea_array; dep_array; cap_array];
            
     abmresults(i).balancesheet = array2table(BS_array);
     abmresults(i).balancesheet.Properties.RowNames = BS_variables;

% Investment
    Inv_array = [banks(i).balancesheet.assets.des_investment; banks(i).balancesheet.assets.investment];
    
    abmresults(i).investment = array2table(Inv_array);
    abmresults(i).investment.Properties.RowNames = Inv_variables;
    
% Interbank market


    lender_array =   [linspace(1,length(banks(i).IBM.L_tot_requests),length(banks(i).IBM.L_tot_requests));
                                                    banks(i).IBM.NB;
                                                    banks(i).IBM.L_tot_requests;
                                                    (banks(i).balancesheet.assets.cash(1:banks(i).failtime,2).*(1-MRR))';
                                                    banks(i).IBM.L_prov_tot_loans;
                                                    banks(i).IBM.hoarding;
                                                    banks(i).IBM.L_tot_loans;
                                                    banks(i).IBM.L_tot_exp_repay;
                                                    banks(i).IBM.L_tot_repaid_loans;
                                                    banks(i).balancesheet.assets.cash(1:banks(i).failtime,3)'];
                                                
   tempL =  lender_array;
   clearvars lender_array
    
   lender_array = tempL(:,all(~isnan(tempL)));
   clearvars tempL
   
   if ~isempty(lender_array)
       abmresults(i).interbankmarket.lenderactivity = array2table(lender_array(2:end,:));
       abmresults(i).interbankmarket.lenderactivity.Properties.RowNames      = LA_variables;
       abmresults(i).interbankmarket.lenderactivity.Properties.VariableNames = strcat('t',strtrim(cellstr(num2str(lender_array(1,:)'))'));
   end
                                              
   borrower_array = [linspace(1,length(banks(i).IBM.B_tot_requests),length(banks(i).IBM.B_tot_requests))
                                                    banks(i).IBM.NL;
                                                    banks(i).IBM.B_tot_requests;
                                                    banks(i).IBM.B_tot_loans;
                                                    banks(i).IBM.B_req_tot_loanrepay;
                                                    (banks(i).balancesheet.assets.cash(1:banks(i).failtime,2).*(1-MRR))';
                                                    banks(i).firesales.tot_des_FS;
                                                    banks(i).firesales.final_firesales;
                                                    banks(i).IBM.B_fin_tot_loanrepay
                                                    banks(i).balancesheet.assets.cash(1:banks(i).failtime,3)'];
                                              
   tempB = borrower_array;
   clearvars borrower_array
    
   borrower_array = tempB(:,all(~isnan(tempB)));
   clearvars tempB
   
   if ~isempty(borrower_array)
        abmresults(i).interbankmarket.borroweractivity = array2table(borrower_array(2:end,:));
        abmresults(i).interbankmarket.borroweractivity.Properties.RowNames      = BA_variables;
        abmresults(i).interbankmarket.borroweractivity.Properties.VariableNames = strcat('t',strtrim(cellstr(num2str(borrower_array(1,:)'))'));
   end
   
   abmresults(i).failtime = banks(i).failtime;
 
end
end
