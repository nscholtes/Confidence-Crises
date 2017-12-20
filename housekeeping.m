function [banks] = housekeeping(banks,n_banks)

 if (isfield(banks,'lending_cps') && isfield(banks,'borrowing_cps'))
    
        remove_IB_vars_1 = {'lending_cps', 'nonlending_cps', 'borrowing_cps','nonborrowing_cps'};
    
        banks = rmfield(banks,remove_IB_vars_1);
 end
    
remove_IB_vars_2 = 'counterpartyids';
    
banks = rmfield(banks,remove_IB_vars_2);

remove_IB_vars_3 = {'B_bil_requests', 'B_bil_loans', 'B_req_bil_loanrepay', 'B_fin_bil_loanrepay','L_bil_loans',...
    'L_bil_requests','L_bil_exp_repay','L_bil_repaid_loans'};

remove_IB_vars_4 = {'des_FS_vec','act_FS_vec'};

for i = 1:n_banks
    for j = 1:numel(remove_IB_vars_3)
            banks(i).IBM.(remove_IB_vars_3{j}) = [];
    end
    
    for j = 1:numel(remove_IB_vars_4)
            banks(i).firesales.(remove_IB_vars_4{j}) = [];
    end
end
    
end

