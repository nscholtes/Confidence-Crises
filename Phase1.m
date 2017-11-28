function [banks,DBV,DLV,NT_matrices,NNT_matrices]...
    =Phase1(banks,ActiveBanks,ibn_adjmat,n_banks,beta,Pars_interestrates,Pars_dshock,MRR,...
    NT_matrices,NNT_matrices,t,fileID_D,writeoption,tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Phase 1: Deposit shock and move to interbank markets. Conditional on the
% SIGN of the deposit shock, banks are classified either as lenders 
% or borrowers on the interbank market and activate their lender or borrower
% connections as a result of this.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------
%% INITIALISATION
%-------------------------------------------------------------------------

tau = 2;

% Assigning variables

r_d  = Pars_interestrates(2);
r_e  = Pars_interestrates(3);
r_b  = Pars_interestrates(4);

theta = Pars_dshock(1);
SPF   = Pars_dshock(2); 

% Initialising counters used in for-loops

lender_count = 0;
borr_count   = 0;

B_noIB_count = 0;
B_noIB       = [];

% Initialising vectors and matrices 

% NxT matrices

shocksign_mat         = NT_matrices(:,1:t,1);
designated_lender_mat = NT_matrices(:,t,2);
designated_borr_mat   = NT_matrices(:,t,3);

% NxNxT matrices

expost_adjmat         = NNT_matrices(:,:,t,1);
borr_req_mat          = NNT_matrices(:,:,1:t,2);
borr_rec_mat          = NNT_matrices(:,:,t,3);
borr_rep_mat          = NNT_matrices(:,:,1:t,4);
lender_lend_mat       = NNT_matrices(:,:,t,5);

adj_lender_lend_mat   = NNT_matrices(:,:,1:t,6);
RP_mat                = NNT_matrices(:,:,t,7);
IBrate_mat            = NNT_matrices(:,:,1:t,8);

%-------------------------------------------------------------------------------------
%% DEPOSIT SHOCK LOOP
%-------------------------------------------------------------------------------------   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametrise the deposit shock based on past bank deposits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dshock_store = zeros(1,numel(ActiveBanks));

for i = ActiveBanks
    
    epsilon_t = normrnd(0,1);
    banks(i).depositshock(t) = theta*((beta*banks(i).balancesheet.assets.total(t,tau-1))-banks(i).balancesheet.liabilities.deposits(t,tau-1))+...
        SPF*banks(i).balancesheet.liabilities.deposits(t,tau-1)*epsilon_t;
    
    dshock_store(i) =  banks(i).depositshock(t);
            
end

dshock_mean = mean(dshock_store);

for i = ActiveBanks
    % Update deposits
    
    banks(i).depositshock(t) = banks(i).depositshock(t) - dshock_mean;
    
    banks(i).balancesheet.liabilities.deposits(t,tau) = banks(i).balancesheet.liabilities.deposits(t,tau-1) + banks(i).depositshock(t);
     %banks(i).balancesheet.liabilities.deposits(t,tau) = banks(i).depositshock(t);
    
    banks(i).balancesheet.assets.cash(t,tau) = banks(i).balancesheet.assets.cash(t,tau-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assigning to banks facing a POSITIVE shock the status of 'designated lender'
% --> These banks redirect the additional liquidity to finance an illiquid investment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
    if  banks(i).depositshock(t) >= 0
           
        banks(i).IBM.status(t) = 'L';
        banks(i).IBM.statusvec(t) = 1;
        shocksign_mat(i,t) = 1;
        
         banks(i).balancesheet.assets.cash(t,tau) = banks(i).balancesheet.assets.cash(t,tau)...
         + r_e.* banks(i).balancesheet.assets.external_assets(t,tau-1)... % Receive interest on external asset portfolio;
         - r_d.*banks(i).balancesheet.liabilities.deposits(t,tau); % Interest on new deposits paid out
        
        if t == 1
            banks(i).balancesheet.assets.des_investment(t) = banks(i).depositshock(t);
            banks(i).balancesheet.assets.investment(t)     = banks(i).balancesheet.assets.des_investment(t);
            
            %banks(i).balancesheet.assets.cash(t,tau)       = banks(i).balancesheet.assets.cash(t,tau-1);
            
        elseif t>1
            
            banks(i).balancesheet.assets.des_investment(t)     =  nanmean(banks(i).balancesheet.assets.investment(1:t-1));
            
% Banks seek to maintain a constant level of private sector investment in each period.

% CASE I: Deposit shock sufficient to cover investment requirement (best-case scenario)
            if banks(i).depositshock(t) >= banks(i).balancesheet.assets.des_investment(t)
                
                banks(i).balancesheet.assets.investment(t) = banks(i).balancesheet.assets.des_investment(t);
                
                % 1st end-of-phase update (LENDERS): Add remaining deposits to reserves after investment
                banks(i).balancesheet.assets.cash(t,tau) = banks(i).balancesheet.assets.cash(t,tau) ...
                    + (banks(i).depositshock(t)  - banks(i).balancesheet.assets.investment(t)); 
                                
% CASE II: Deposit shock insufficient but difference can be covered by excess reserves (intermediate scenario) 
%%% Make investment but results in a decrease in reserves
            elseif banks(i).depositshock(t) < banks(i).balancesheet.assets.des_investment(t)  &&...
                    (banks(i).balancesheet.assets.des_investment(t)  - banks(i).depositshock(t))  <= (1-MRR)*banks(i).balancesheet.assets.cash(t,tau)
                
                 banks(i).balancesheet.assets.investment(t)     =  banks(i).balancesheet.assets.des_investment(t);
                 
                 banks(i).balancesheet.assets.cash(t,tau) = banks(i).balancesheet.assets.cash(t,tau)...
                     + (banks(i).depositshock(t) - banks(i).balancesheet.assets.investment(t));
                                  
  % CASE III: Deposit shock insufficient to make desired investment and excess reserves insufficient to cover difference (worst-case scenario)  
  %%% Allocate up to maximum allowable (such that minimum reserve requirement becomes binding)
             elseif banks(i).depositshock(t) < banks(i).balancesheet.assets.des_investment(t)  &&...
                    (banks(i).balancesheet.assets.des_investment(t)  - banks(i).depositshock(t))  > (1-MRR)*banks(i).balancesheet.assets.cash(t,tau)
                
                 banks(i).balancesheet.assets.investment(t)     =  (1-MRR)*banks(i).balancesheet.assets.cash(t,tau)+banks(i).depositshock(t);
                 
                 banks(i).balancesheet.assets.cash(t,tau) = (MRR)*banks(i).balancesheet.assets.cash(t,tau);
            else
                disp('Error in lender investment decision!')
                                     
            end
                    
% Populating matrix with all designated borrowers in each period

         end
        
        lender_count = lender_count+1;
        designated_lender_mat(i) = banks(i).id;
      
% Creating a matrix with all designated lenders in each period

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assigning to banks facing a NEGATIVE shock the status of 'designated
% borrower' --> These banks, wishing to make an investment, must overcome
% the negative shock, which is financed by INTERBANK BORROWING --> TOTAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif banks(i).depositshock(t) < 0

        banks(i).IBM.status(t)    = 'B';
        banks(i).IBM.statusvec(t) = -1;
        shocksign_mat(i,t)        = -1;
    
        borr_count             = borr_count+1;
        designated_borr_mat(i) = banks(i).id; 
        
    end
end
 
%-------------------------------------------------------------------------
%% USING SHOCK SIGN TO UPDATE ADJACENCY MATRIX
%-------------------------------------------------------------------------

temp_DL = designated_lender_mat;
temp_DB = designated_borr_mat;

% Initialising final unweighted, directed adjacency matrix containing the interbank
% market transactions occuring in that period due to shock sign
% distribution

expost_adjmat(:,:) = ibn_adjmat;

% Removing zeros from lender/borrower matrix --> Vectors containing the ids
% of designated lenders (DLV) and designated borrowers (DBV)

DLV = temp_DL(temp_DL ~=0);
DBV = temp_DB(temp_DB ~=0);

% Loop through the vector of designated borrowers to get rid of other borrowers in its neighbourhood
%%% Applied when requesting loans (depending on calibrated 'interbank
%%% market transparency')

for i = 1:numel(DBV)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Defines the vector of 'eligible lending counterparties' for each borrower i.e. banks
% that were hit by a positive shock in the current period and can actually lend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elig_cps   = banks(DBV(i)).counterpartyids;
    inelig_cps = intersect(banks(DBV(i)).counterpartyids,DBV)';
    overlap    = ismember(elig_cps, inelig_cps);
    index_inel = find(overlap);
    index_el   = find(~overlap);
    elig_cps(index_inel) = [];
    banks(DBV(i)).lending_cps = elig_cps;
    banks(DBV(i)).nonlending_cps = inelig_cps;
    banks(DBV(i)).numnonlending_cps(t) =  banks(DBV(i)).num_counterparties(t)-numel(banks(DBV(i)).lending_cps);
    
    banks(DBV(i)).numlending_cps(t) = numel(banks(DBV(i)).lending_cps);
    
    banks(DBV(i)).IBM.NL(t) = banks(DBV(i)).numlending_cps(t);
    banks(DBV(i)).IBM.NB(t) = NaN;
    
    %-/+ Filter: Borrowers do not lend to lenders
    
    expost_adjmat(DBV(i),banks(DBV(i)).counterpartyids) = 0;  
    
    % -/- Filter: borrowers will not request liquidity from other borrowers under perfect information:
        
    expost_adjmat(DBV(i),banks(DBV(i)).nonlending_cps) = zeros(1,numel(banks(DBV(i)).nonlending_cps));

end

for i = 1:numel(DLV)
    
    elig_cps   = banks(DLV(i)).counterpartyids;
    inelig_cps = intersect(banks(DLV(i)).counterpartyids,DLV)'; % Other lenders are ineligible for loans
    overlap    = ismember(elig_cps, inelig_cps);
    index_inel = find(overlap);
    index_el   = find(~overlap);
    elig_cps(index_inel) = [];
    banks(DLV(i)).borrowing_cps = elig_cps;
    banks(DLV(i)).nonborrowing_cps = inelig_cps;
    banks(DLV(i)).numborrowing_cps(t) = numel(banks(DLV(i)).borrowing_cps);
    
    banks(DLV(i)).IBM.NB(t) = banks(DLV(i)).numborrowing_cps(t);
    banks(DLV(i)).IBM.NL(t) = NaN;
    
    %+/+ Filter: Remove lender-lender edges
     
    expost_adjmat(DLV(i),banks(DLV(i)).nonborrowing_cps) = zeros(1,numel(banks(DLV(i)).nonborrowing_cps)); 
              
end

%-------------------------------------------------------------------------
%% BORROWER REQUESTS
%-------------------------------------------------------------------------

myprint(writeoption,fileID_D,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\r\n');
myprint(writeoption,fileID_D,'--------------------------------------------------------------------- Period %d Output -----------------------------------------------------------------------\r\n',t);
myprint(writeoption,fileID_D,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\r\n');
myprint(writeoption,fileID_D,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');
myprint(writeoption,fileID_D,'Phase 1\r\n');
myprint(writeoption,fileID_D,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\r\n');
myprint(writeoption,fileID_D,'**********************************************************************************\r\n');
%myprint(writeoption,fileID_D,'**************** Information in the current simulation is %s ****************\r\n',information);
myprint(writeoption,fileID_D,'**********************************************************************************\r\n');

myprint(writeoption,fileID_D,'====================================================================================================================================================================\r\n');
myprint(writeoption,fileID_D,'Step 1: Borrower requests\r\n');
myprint(writeoption,fileID_D,'====================================================================================================================================================================\r\n');

for i = 1:numel(DBV)
    
    % Variables associated to lenders are set = NaN
    
    banks(DBV(i)).numborrowing_cps(t)     = NaN;
    
    banks(DBV(i)).IBM.L_noCP(t)           = NaN;
    banks(DBV(i)).IBM.L_bil_requests      = NaN;
    banks(DBV(i)).IBM.L_tot_requests(t)   = NaN;
    
    banks(DBV(i)).IBM.L_bil_loans(t)      = NaN;
    
    banks(DBV(i)).IBM.L_prov_tot_loans(t) = NaN;
    banks(DBV(i)).IBM.L_tot_loans(t)      = NaN;
    
    banks(DBV(i)).IBM.hoarding(t)             = NaN;
    banks(DBV(i)).IBM.L_hoardingmultiplier(t) = NaN;

% Determining  borrowers' desired investment as a function of shock size
% (period 1) or average past investments (periods > 1)

if t == 1
    banks(DBV(i)).balancesheet.assets.des_investment(t) = abs(banks(DBV(i)).depositshock(t));
elseif t > 1
    banks(DBV(i)).balancesheet.assets.des_investment(t) =  nanmean(banks(DBV(i)).balancesheet.assets.investment(1:t-1));  
end

    banks(DBV(i)).balancesheet.assets.cash(t,tau) = banks(DBV(i)).balancesheet.assets.cash(t,tau)...
         + r_e.* banks(DBV(i)).balancesheet.assets.external_assets(t,tau-1)... % Receive interest on external asset portfolio;
         - r_d.*banks(DBV(i)).balancesheet.liabilities.deposits(t,tau);        % Interest on new deposits paid out
              
% Contingency when shock sign distribution is such that borrower i has no counterparties from which it can request liquidity 
    
    if isempty(banks(DBV(i)).lending_cps)
        
        myprint(writeoption,fileID_D,'Borrower %d has NO lending counterparties in period %d\r\n',DBV(i),t);
        myprint(writeoption,fileID_D,'----------------------------------------------------------------------------------\r\n');
        
        banks(DBV(i)).IBM.B_bil_requests = 0;
        banks(DBV(i)).balancesheet.liabilities.IB_borrowing(t) = 0;
        borr_req_mat(DBV(i),:,t) = zeros(1,n_banks);
     
        banks(DBV(i)).IBM.B_tot_requests(t) = 0;
        banks(DBV(i)).IBM.B_tot_loans(t)    = 0;
        banks(DBV(i)).IBM.B_noCP(t)         = 1;
        
        if (1-MRR)*banks(DBV(i)).balancesheet.assets.cash(t,tau) >= banks(DBV(i)).balancesheet.assets.des_investment(t)
            
             banks(DBV(i)).balancesheet.assets.investment(t) = banks(DBV(i)).balancesheet.assets.des_investment(t);
             
            myprint(writeoption,fileID_D,'--> Sufficient reserves to make investment (%.3f>%.3f)\r\n',...
            banks(DBV(i)).balancesheet.assets.cash(t,tau),banks(DBV(i)).balancesheet.assets.des_investment(t));
             
             banks(DBV(i)).balancesheet.assets.cash(t,tau)   = banks(DBV(i)).balancesheet.assets.cash(t,tau) - ...
                 banks(DBV(i)).balancesheet.assets.investment(t);
            
            myprint(writeoption,fileID_D,'--> Updated reserves given by %.3f\r\n',banks(DBV(i)).balancesheet.assets.cash(t,tau));
        
        elseif (1-MRR)*banks(DBV(i)).balancesheet.assets.cash(t,tau) < banks(DBV(i)).balancesheet.assets.des_investment(t)
            
            banks(DBV(i)).balancesheet.assets.investment(t) = (1-MRR)*banks(DBV(i)).balancesheet.assets.cash(t,tau);
            
            myprint(writeoption,fileID_D,'--> Insufficient reserves (%.3f<%3f) to make investment\r\n',...
            (1-MRR)*banks(DBV(i)).balancesheet.assets.cash(t,tau),banks(DBV(i)).balancesheet.assets.des_investment(t));
            
            banks(DBV(i)).balancesheet.assets.cash(t,tau)   = MRR*banks(DBV(i)).balancesheet.assets.cash(t,tau);
            
            myprint(writeoption,fileID_D,'Allocate up to maximum allowable. Updated reserves given by %.3f\r\n',...
                banks(DBV(i)).balancesheet.assets.cash(t,tau));
        else
            disp('Error for borrowers with sufficient cash to make investment')
        end 
          
    elseif ~isempty(banks(DBV(i)).lending_cps)
    
        myprint(writeoption,fileID_D,'Borrower %d has %d lending counterparties in period %d\r\n',DBV(i),banks(DBV(i)).numlending_cps(t),t);
        myprint(writeoption,fileID_D,'--> Borrower %d has %d counterparties who are borrowers themselves in period %d\r\n',DBV(i),banks(DBV(i)).numnonlending_cps(t),t);
        
%          Case I: Can make desire investment using excess reserves --> No need for interbank borrowing
%         if banks(DBV(i)).balancesheet.assets.cash(t,tau)*(1-MRR) >= banks(DBV(i)).balancesheet.assets.des_investment(t)
%             
%             B_noIB_count         = B_noIB_count+1;
%             B_noIB(B_noIB_count) = DBV(i);
%             
%             banks(DBV(i)).IBM.B_tot_requests(t) = 0;
%             banks(DBV(i)).IBM.B_bil_requests    =  zeros(1,banks(DBV(i)).numlending_cps(t));
%             
%             banks(DBV(i)).IBM.B_bil_requests =  ones(1,banks(DBV(i)).numlending_cps(t))...
%             .*((1./banks(DBV(i)).numlending_cps(t)).*banks(DBV(i)).IBM.B_tot_requests(t));
%         
%             borr_req_mat(DBV(i),banks(DBV(i)).lending_cps,t) = banks(DBV(i)).IBM.B_bil_requests;
%         
%             banks(DBV(i)).IBM.B_noCP(t) = 0;
%             
%             banks(DBV(i)).balancesheet.assets.investment(t) = banks(DBV(i)).balancesheet.assets.des_investment(t);
%             
%             banks(DBV(i)).balancesheet.assets.cash(t,tau) =  banks(DBV(i)).balancesheet.assets.cash(t,tau)-...
%                 banks(DBV(i)).balancesheet.assets.investment(t);
%             
%         Case II: Inadequate reserves: Allocate up to maximum allowable and obtain remainder from the interbank market     
%         elseif banks(DBV(i)).balancesheet.assets.cash(t,tau)*(1-MRR) < banks(DBV(i)).balancesheet.assets.des_investment(t)
            
            banks(DBV(i)).IBM.B_tot_requests(t) = banks(DBV(i)).balancesheet.assets.des_investment(t) ;
            
            %-...
               % banks(DBV(i)).balancesheet.assets.cash(t,tau)*(1-MRR) ;
        
            banks(DBV(i)).IBM.B_bil_requests = ones(1,banks(DBV(i)).numlending_cps(t))...
            .*((1./banks(DBV(i)).numlending_cps(t)).*banks(DBV(i)).IBM.B_tot_requests(t));
        
            %banks(DBV(i)).balancesheet.assets.cash(t,tau) =  MRR*banks(DBV(i)).balancesheet.assets.cash(t,tau);
        
            borr_req_mat(DBV(i),banks(DBV(i)).lending_cps,t) = banks(DBV(i)).IBM.B_bil_requests;
        
            banks(DBV(i)).IBM.B_noCP(t) = 0;
        
            myprint(writeoption,fileID_D,'Borrower %d requests %.3f from each  of the %d lending counterparties in period %d\r\n',...
                DBV(i),(1./banks(DBV(i)).numlending_cps(t)).*banks(DBV(i)).IBM.B_tot_requests(t),banks(DBV(i)).numlending_cps(t),t);
            myprint(writeoption,fileID_D,'----------------------------------------------------------------------------------\r\n');
        
        %end
    end
end
% For consistency: transpose such that rows (columns) correspond to lender (borrowers)

borr_req_mat(:,:,t) = borr_req_mat(:,:,t)';

% Remove borrowers able to cover liquidity needs through reallocation of reserves
% --> Do not go to interbank market for loans

% if B_noIB_count > 0
%     [del_no_IB, ~, ~] = find(DBV==B_noIB);
%     DBV(del_no_IB)    = [];
%     
%     for i = B_noIB
%     
%     banks(i).IBM.L_bil_exp_repay    = NaN;
%     banks(i).IBM.L_tot_exp_repay(t) = NaN;
%     
%     banks(i).IBM.B_bil_loans    = zeros(1,banks(i).numlending_cps(t));
%     banks(i).IBM.B_tot_loans(t) = 0;
%                     
%     end
%         
%     for i=1:numel(DLV)
%             
%         temp_borrowers = banks(DLV(i)).borrowing_cps
%         
%         [del_TB, ~, ~] = find(temp_borrowers' == B_noIB);
%         temp_borrowers(del_TB)    = [];
%         
%         temp_borrowers
%             
%         banks(DLV(i)).borrowing_cps       = temp_borrowers;
%         banks(DLV(i)).numborrowing_cps(t) = numel(banks(DLV(i)).borrowing_cps);
%             
%         %if isempty(banks(DLV(i)).borrowing_cps)
%           % L_noIB = DLV(i);           
%        %end   
%     
%     end  
%     %[del_DLV, ~, ~] = find(DLV==L_noIB);
%     %DLV(del_DLV)    = [];
% end


%-------------------------------------------------------------------------
%% LENDER OFFERS
%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lenders redistribute available liquidity equally across borrowing neighbours
%%% Banks will meet the aggregate request provided  that their reserves do
%%% not fall below the buffer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Redistribution of available liquidity (cash reserves+interest on external
% assets - interest payments to depositors) across borrowing counterparties

myprint(writeoption,fileID_D,'====================================================================================================================================================================\r\n');
myprint(writeoption,fileID_D,'Step 2: Lender offers\r\n');
myprint(writeoption,fileID_D,'====================================================================================================================================================================\r\n');

for i = 1:numel(DLV)
    
    % Variables associated to borrowers are set = NaN
    
    banks(DLV(i)).numlending_cps(t)     = NaN;
    
    banks(DLV(i)).IBM.B_noCP(t)         = NaN;
    banks(DLV(i)).IBM.B_suff_loans(t)   = NaN;
    
    banks(DLV(i)).IBM.B_bil_requests    = NaN;
    banks(DLV(i)).IBM.B_tot_requests(t) = NaN;
    
    banks(DLV(i)).IBM.B_bil_loans(t)    = NaN;
    banks(DLV(i)).IBM.B_tot_loans(t)    = NaN;   
    
    % 2nd end of phase update (LENDERS): Lenders determine cash available + extra income to make interbank loans
      
    temp_reqstore = borr_req_mat(DLV(i),:,t);
    temp_reqstore(temp_reqstore==0)=[];
    
    banks(DLV(i)).IBM.L_bil_requests    = temp_reqstore;    
    banks(DLV(i)).IBM.L_tot_requests(t) = sum(banks(DLV(i)).IBM.L_bil_requests);
    
    banks(DLV(i)).borrowing_cps       = setdiff(banks(DLV(i)).borrowing_cps,B_noIB);
    banks(DLV(i)).numborrowing_cps(t) = numel(banks(DLV(i)).borrowing_cps);
        
    clearvars temp_reqstore
   
    if isempty(banks(DLV(i)).borrowing_cps)
        
        myprint(writeoption,fileID_D,'Lender %d has NO borrowing counterparties in period %d\r\n',DLV(i),t);
        myprint(writeoption,fileID_D,'----------------------------------------------------------------------------------\r\n');
          
        banks(DLV(i)).IBM.hoarding(t)             = 0; % Lender hoarding
        banks(DLV(i)).IBM.L_hoardingmultiplier(t) = 0;
        banks(DLV(i)).IBM.L_bil_loans             = 0;
        banks(DLV(i)).IBM.L_prov_tot_loans(t)     = 0;
        banks(DLV(i)).IBM.L_tot_loans(t)          = 0;
        
        banks(DLV(i)).IBM.L_noCP(t) = 1;
        
        lender_lend_mat(DLV(i),:) = zeros(1,n_banks);
        
    else
        
        myprint(writeoption,fileID_D,'Lender %d has %d counterparties\r\n',DLV(i),banks(DLV(i)).num_counterparties(t));
        myprint(writeoption,fileID_D,'--> Lender %d has %d borrowing counterparties in period %d\r\n',DLV(i),banks(DLV(i)).numborrowing_cps(t),t);
        myprint(writeoption,fileID_D,'Lender %d received loan requests =  %.3f  in period %d\r\n',DLV(i),banks(DLV(i)).IBM.L_tot_requests(t),t);
        
        banks(DLV(i)).IBM.L_noCP(t)               = 0;
        banks(DLV(i)).IBM.hoarding(t)             = 0;
        banks(DLV(i)).IBM.L_hoardingmultiplier(t) = 0;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total lending step 1: Compare total requests with available cash reserves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % prov_tot_loans: loan amounts prior to hoarding (but has to satisfy banks' available liquidity)
        banks(DLV(i)).IBM.L_prov_tot_loans(t) = min(banks(DLV(i)).balancesheet.assets.cash(t,tau)*(1-MRR),banks(DLV(i)).IBM.L_tot_requests(t));
        
        myprint(writeoption,fileID_D,'*1.CHECK AVAILABLE LIQUIDITY*\r\n');
        
        % CASE I: Total requests <= Available reserves --> Lenders makes all requested loans.
            
            if abs(banks(DLV(i)).IBM.L_prov_tot_loans(t) - banks(DLV(i)).IBM.L_tot_requests(t)) < tol
                
                banks(DLV(i)).IBM.L_tot_loans(t) = banks(DLV(i)).IBM.L_prov_tot_loans(t);
                                
                banks(DLV(i)).IBM.L_bil_loans = banks(DLV(i)).IBM.L_bil_requests;
                           
                myprint(writeoption,fileID_D,'Lender %d has sufficient reserves (%.3f) to match all requests (%.3f) in period %d\r\n',DLV(i),...
                    banks(DLV(i)).balancesheet.assets.cash(t,tau),banks(DLV(i)).IBM.L_tot_requests(t),t);    
                myprint(writeoption,fileID_D,'----------------------------------------------------------------------------------\r\n');
                
         % CASE II:  Total requests > available reserves --> Lender allocates maximum allowable which is divided equally across counterparties and
         % compared to each bilateral request individually using the same min{allocated bilateral funds;request} approach.
               
            elseif abs(banks(DLV(i)).IBM.L_prov_tot_loans(t) - banks(DLV(i)).balancesheet.assets.cash(t,tau)*(1-MRR)) < tol
                
                %banks(DLV(i)).IBM.L_tot_loans(t) == banks(DLV(i)).balancesheet.assets.cash(t,tau)*(1-MRR)
                                
                    for j = 1:banks(DLV(i)).numborrowing_cps(t)
                    
                        banks(DLV(i)).IBM.L_bil_loans(j) = min(banks(DLV(i)).IBM.L_bil_requests(j),...
                            (banks(DLV(i)).balancesheet.assets.cash(t,tau)*(1-MRR))/(banks(DLV(i)).numborrowing_cps(t)));
                        
                    end
                        
                        banks(DLV(i)).IBM.L_tot_loans(t) = sum(banks(DLV(i)).IBM.L_bil_loans);
                        
                        myprint(writeoption,fileID_D,'Total requests (%.3f) to Lender %d exceed available reserves (%.3f). Allocate %.3f in period %d\r\n',...
                            banks(DLV(i)).IBM.L_tot_requests(t),DLV(i),...
                            banks(DLV(i)).balancesheet.assets.cash(t,tau)*(1-MRR),banks(DLV(i)).IBM.L_tot_loans(t),t);
                        myprint(writeoption,fileID_D,'----------------------------------------------------------------------------------\r\n');  
            else
                disp('Error in lender reserve allocation')
            end  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Total lending step 2: Use past funding conditions to determine current level of loans
 %%% Search through interbank history to last period where I was liquidity short and had to request lending. Compare total
 %%% requests to total received loans ==> Proxy for banks' perception of funding liquidity risk
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        if t>1      
            t_search = t-1;
                while t_search >= 1
                    
                    if shocksign_mat(DLV(i),t_search) == -1 && banks(DLV(i)).IBM.B_tot_requests(t_search)~=0 ...
                            && banks(DLV(i)).IBM.B_tot_loans(t_search)~=0 && banks(DLV(i)).IBM.B_tot_requests(t_search)~=0
                        
                        myprint(writeoption,fileID_D,'*2. HOARDING DECISION*\r\n');
                        
                        banks(DLV(i)).IBM.L_hoardingmultiplier(t) = (banks(DLV(i)).IBM.B_tot_loans(t_search)./...
                            banks(DLV(i)).IBM.B_tot_requests(t_search));
                                         
                        banks(DLV(i)).IBM.L_tot_loans(t)  = banks(DLV(i)).IBM.L_hoardingmultiplier(t).*banks(DLV(i)).IBM.L_tot_loans(t);
                        
                        banks(DLV(i)).IBM.hoarding(t)   = (1-banks(DLV(i)).IBM.L_hoardingmultiplier(t)).*banks(DLV(i)).IBM.L_tot_loans(t);
                        
                        if abs(banks(DLV(i)).IBM.hoarding(t)) < tol
                            banks(DLV(i)).IBM.hoarding(t) = 0;
                        end
                        
                        myprint(writeoption,fileID_D,'Lender %d was a borrower in period %d, receiving %.3f percent of total requested loans\r\n',...
                            DLV(i),t_search,banks(DLV(i)).IBM.B_tot_loans(t_search)/banks(DLV(i)).IBM.B_tot_requests(t_search)*100);
                        
                        % Case I: Total loans = total requests after hoarding => Each individial request is met
                        if abs(banks(DLV(i)).IBM.L_tot_loans(t) - banks(DLV(i)).IBM.L_tot_requests(t))< tol
                            
                             %banks(DLV(i)).IBM.L_hoardingmultiplier(t).*
                            %banks(DLV(i)).IBM.L_tot_loans(t) == banks(DLV(i)).IBM.L_tot_requests(t)
                                            
                            banks(DLV(i)).IBM.L_bil_loans = banks(DLV(i)).IBM.L_bil_requests;
                            
                            banks(DLV(i)).IBM.L_tot_loans(t) = sum(banks(DLV(i)).IBM.L_bil_loans);

                            myprint(writeoption,fileID_D,'--> After hoarding, Lender %d has sufficient reserves to match all requests in period %d\r\n',DLV(i),t);
                            
                        % Case II: Total requests > total loans after hoarding => Allocated total loan amount allocated on a
                        % case by case basis.
                        elseif  abs(banks(DLV(i)).IBM.L_tot_loans(t) - banks(DLV(i)).IBM.L_tot_requests(t)) >= tol
                            
                            %if abs(banks(DLV(i)).IBM.L_tot_loans(t) -  banks(DLV(i)).IBM.L_hoardingmultiplier(t).*banks(DLV(i)).balancesheet.assets.cash(t,tau)*(1-MRR)) < tol
                            
                           % banks(DLV(i)).IBM.L_tot_loans(t) == banks(DLV(i)).balancesheet.assets.cash(t,tau)*(1-MRR)
                                                    
                            for j = 1:banks(DLV(i)).numborrowing_cps(t)
                    
                                banks(DLV(i)).IBM.L_bil_loans(j) = min(banks(DLV(i)).IBM.L_bil_requests(j),...
                                    (banks(DLV(i)).IBM.L_tot_loans(t))/(banks(DLV(i)).numborrowing_cps(t)));
                            end
                            
                            myprint(writeoption,fileID_D,'--> After hoarding, Total requests to Lender %d exceed available reserves. Allocate %.3f in period %d\r\n',...
                                DLV(i),banks(DLV(i)).IBM.L_tot_loans(t),t);
                            
                             banks(DLV(i)).IBM.L_tot_loans(t) = sum(banks(DLV(i)).IBM.L_bil_loans);
                             
                        else
                            disp('Error in bank hoarding decision!')
                             
                        end
                             
                        myprint(writeoption,fileID_D,'Lender %d lent out %.3f and hoarded %.3f in period %d\r\n',...
                            DLV(i),banks(DLV(i)).IBM.L_tot_loans(t),banks(DLV(i)).IBM.hoarding(t),t);
                        
                        for j = 1:banks(DLV(i)).numborrowing_cps(t) 
                            myprint(writeoption,fileID_D,'--> Borrower %d requested %.3f and received %.3f from Lender %d in period %d\r\n' ,...
                            banks(DLV(i)).borrowing_cps(j),banks(DLV(i)).IBM.L_bil_requests(j),banks(DLV(i)).IBM.L_bil_loans(j),DLV(i),t);
                        end   

                        t_search = 0;   % For now, only look at most recent period where lender was liquidity short
                    end                  
                    t_search = t_search -1; 
                                
                end
                myprint(writeoption,fileID_D,'----------------------------------------------------------------------------------\r\n');
        end
        
    lender_lend_mat(DLV(i),banks(DLV(i)).borrowing_cps) = banks(DLV(i)).IBM.L_bil_loans;

    end
                         
    % 3rd end of phase update (LENDERS): Subtract total interbank loans from reserves
        
    banks(DLV(i)).balancesheet.assets.cash(t,tau) = banks(DLV(i)).balancesheet.assets.cash(t,tau)-banks(DLV(i)).IBM.L_tot_loans(t);

end

borr_rec_mat = lender_lend_mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transferring lender offers to borrower classes
%%% Borrowers observe how much of the requested liquidity they were able to obtain on the interbank market
%%% Use liquidity to make investment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myprint(writeoption,fileID_D,'====================================================================================================================================================================\r\n');
myprint(writeoption,fileID_D,'Step 3: Borrower response to offers\r\n');
myprint(writeoption,fileID_D,'====================================================================================================================================================================\r\n');

for i = 1:numel(DBV)
    
    banks(DBV(i)).IBM.L_bil_exp_repay    = NaN;
    banks(DBV(i)).IBM.L_tot_exp_repay(t) = NaN;
    
    % 2nd end of phase update (BORROWERS): Borrowers determine available reserves prior to receiving interbank loans
                            
    temp_loanstore = lender_lend_mat(:,DBV(i))';
    temp_loanstore(temp_loanstore==0)=[];
    
    banks(DBV(i)).IBM.B_bil_loans = temp_loanstore;
    
    clearvars temp_loanstore
  
    banks(DBV(i)).IBM.B_tot_loans(t)   = sum(banks(DBV(i)).IBM.B_bil_loans);
      
    if banks(DBV(i)).IBM.B_tot_requests(t) > tol && ~isempty(banks(DBV(i)).lending_cps)
                 
        banks(DBV(i)).balancesheet.assets.investment(t) = banks(DBV(i)).IBM.B_tot_loans(t);
        
        if banks(DBV(i)).balancesheet.assets.investment(t) < banks(DBV(i)).balancesheet.assets.des_investment(t)
            
            if (banks(DBV(i)).balancesheet.assets.des_investment(t) - banks(DBV(i)).balancesheet.assets.investment(t)) <...
                (1-MRR)*banks(DBV(i)).balancesheet.assets.cash(t,tau)
                    banks(DBV(i)).balancesheet.assets.investment(t) = banks(DBV(i)).balancesheet.assets.des_investment(t);
                    banks(DBV(i)).balancesheet.assets.cash(t,tau)   = banks(DBV(i)).balancesheet.assets.cash(t,tau) - ...
                        (banks(DBV(i)).balancesheet.assets.des_investment(t) - banks(DBV(i)).balancesheet.assets.investment(t));
            else
                banks(DBV(i)).balancesheet.assets.investment(t) = (1-MRR)*banks(DBV(i)).balancesheet.assets.cash(t,tau);
                banks(DBV(i)).balancesheet.assets.cash(t,tau)   = (MRR)*banks(DBV(i)).balancesheet.assets.cash(t,tau);
            end
        end
        
        %...+ banks(DBV(i)).balancesheet.assets.cash(t,tau)*(1-MRR)
     
        myprint(writeoption,fileID_D,'Borrower %d obtained the requested loans and made the desired investment of %.3f in period %d\r\n',...
            DBV(i),banks(DBV(i)).balancesheet.assets.investment(t),t);
        
        myprint(writeoption,fileID_D,'Borrower %d loan/request = %.3f percent in period %d\r\n',...
            DBV(i),banks(DBV(i)).IBM.B_tot_loans(t)/banks(DBV(i)).IBM.B_tot_requests(t)*100,t);
        myprint(writeoption,fileID_D,'----------------------------------------------------------------------------------\r\n');
    end
end

%-------------------------------------------------------------------------
%% INTEREST RATE SETTING
%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lenders set risk premium on interbank lending volume decided in STEP 1 to 
%counterparties hit by negative shock. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

myprint(writeoption,fileID_D,'====================================================================================================================================================================\r\n');
myprint(writeoption,fileID_D,'Step 4: Lender interest rate setting\r\n');
myprint(writeoption,fileID_D,'====================================================================================================================================================================\r\n');

% Initialisation: interbank rate set at (calibrated) risk-free rate

for i = 1:numel(DLV)
       
    if  isempty(banks(DLV(i)).borrowing_cps)
                
        banks(DLV(i)).IBM.L_bil_exp_repay     = 0;      
        banks(DLV(i)).IBM.L_tot_exp_repay(t)  = 0;
        
        RP_mat(DLV(i),:)                = zeros(1,n_banks);
        IBrate_mat(DLV(i),:,t)          = zeros(1,n_banks);
        adj_lender_lend_mat(DLV(i),:,t) = zeros(1,n_banks);

        %IBrate_mat(DLV(i),:,t)          = (1+r_b)*expost_adjmat(DLV(i),:);
        %adj_lender_lend_mat(DLV(i),:,t) = expost_adjmat(DLV(i),:);
                
        myprint(writeoption,fileID_D,'Lender %d received no requests in period %d\r\n',DLV(i),t);
        myprint(writeoption,fileID_D,'-------------------------------------------------------------------------------------------\r\n');
    
    elseif ~isempty(banks(DLV(i)).borrowing_cps)       
        
        if t==1
            RP_mat(DLV(i),:)                = expost_adjmat(DLV(i),:);
            
            
            IBrate_mat(DLV(i),:,t)          = (1+r_b)*expost_adjmat(DLV(i),:);
            adj_lender_lend_mat(DLV(i),:,t) = IBrate_mat(DLV(i),:,t).*lender_lend_mat(DLV(i),:);
            
            banks(DLV(i)).IBM.L_bil_exp_repay    = adj_lender_lend_mat(DLV(i),banks(DLV(i)).borrowing_cps,t);
            banks(DLV(i)).IBM.L_tot_exp_repay(t) = sum(adj_lender_lend_mat(DLV(i),banks(DLV(i)).borrowing_cps,t));
            
            for j=1:banks(DLV(i)).numborrowing_cps(t)
                myprint(writeoption,fileID_D,'--> Lender %d sets a counterparty risk premium of %.3f on borrower %d in period %d\r\n',...
                    DLV(i),RP_mat(DLV(i),banks(DLV(i)).borrowing_cps(j))-1,banks(DLV(i)).borrowing_cps(j),t);
                    
                myprint(writeoption,fileID_D,'--> The interbank rate for this transaction is %.2f percent\r\n',...
                    (IBrate_mat(DLV(i),banks(DLV(i)).borrowing_cps(j),t)-1)*100);
            
            end
               myprint(writeoption,fileID_D,'-------------------------------------------------------------------------------------------\r\n');

        end
    
        pastLBcount        = 0; pastLB_Bstore       = []; pastLB_tstore    = []; 
        pastLB_loanstore   = []; pastLB_repstore    = []; pastLB_ratestore = [];
        pastLB_loanstoreav = []; pastLB_repstoreav  = [];
        
        nopastLB_Bstore = [];
        
        CD_pastLB_Bstore = []; freqCD  = []; CD_ids = []; CD_index = [];        
        noCD_ids = []; noCD_index = [];
           
        if t>1
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Counterparty risk premium: Locate last time period where current lender acted as a lender
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Looping through each of lender i's CURRENT borrowers

            for j=1:banks(DLV(i)).numborrowing_cps(t)
                
                t_search = t-1;

                while (t_search >= 1)

                    % STOP Condition: lender i and borrower j were in a lender-borrower relationship in the past
                
                    if shocksign_mat(DLV(i),t_search) == 1 && shocksign_mat(banks(DLV(i)).borrowing_cps(j),t_search) == -1 &&...
                        adj_lender_lend_mat(DLV(i),banks(DLV(i)).borrowing_cps(j),t_search) ~= 0 &&...
                        borr_rep_mat(DLV(i),banks(DLV(i)).borrowing_cps(j),t_search)~=0 
                        
                        pastLBcount = pastLBcount+1;         % Count when STOP condition is met
         
                        pastLB_Bstore(pastLBcount) = banks(DLV(i)).borrowing_cps(j);  % Store borrowing counterparty for which STOP condition above is satisfied
                        pastLB_tstore(pastLBcount) = t_search;                       %  Store period in which STOP condition is satisfied
                                           
                        pastLB_loanstore(pastLBcount) = adj_lender_lend_mat(DLV(i),banks(DLV(i)).borrowing_cps(j),t_search); % Loan|STOP 
                        pastLB_repstore(pastLBcount)  = borr_rep_mat(DLV(i),banks(DLV(i)).borrowing_cps(j),t_search);        % Repayment|STOP
                        pastLB_ratestore(pastLBcount) = IBrate_mat(DLV(i),banks(DLV(i)).borrowing_cps(j),t_search);          % IBrate|STOP 
    
                    end                    
                    t_search = t_search - 1;      
                end                
            end
                        
            if pastLBcount ~= 0 % Case when at least one past L-B relationsbip exists betwen current i and all j's
                
                % Compute averages of loan provision and repayment for use in risk premium calculation
                
                for k = 1:pastLBcount
                    temp_loanstoreav = adj_lender_lend_mat(DLV(i),:,pastLB_tstore(k));
                    temp_repstoreav  = borr_rep_mat(DLV(i),:,pastLB_tstore(k));
                     
                    pastLB_loanstoreav(k) = sum(temp_loanstoreav)./sum(temp_loanstoreav~=0);
                    pastLB_repstoreav(k)  = sum(temp_repstoreav )./sum(temp_repstoreav~=0);                   
                end
                
                %pastLB_loanstoreav
                %pastLB_repstoreav 
                
                %CD = Check Duplicates (multiple occurences of past L-B relationships for the same borrowing counterparty)
                
                CD_pastLB_Bstore = unique(pastLB_Bstore);
                
                freqCD = histc(pastLB_Bstore,CD_pastLB_Bstore); % Intermediary step: Used for next line
                
                CD_ids   = CD_pastLB_Bstore(freqCD>1);  % Identity of duplicate values                
                noCD_ids = CD_pastLB_Bstore(freqCD==1); % Identity of non-duplicate values (only one past L/B relationship)
                
                if ~isequal(CD_pastLB_Bstore,pastLB_Bstore) % Duplicate values present

                    for j = 1:length(CD_ids) % Cycle through all borrowers with multiple past L/B relationships with current lender
                            
                        CD_index(j,:) = find(pastLB_Bstore == CD_ids(j)); % Indices of duplicate values
                        
                        RP_mat(DLV(i),CD_ids(j)) = mean(ones(1,numel(pastLB_tstore(CD_index(j,:))))...
                            +(pastLB_loanstore(CD_index(j,:))-pastLB_repstore(CD_index(j,:)))...
                            .*(pastLB_loanstoreav(CD_index(j,:))/pastLB_repstoreav(CD_index(j,:))));
                        
                         IBrate_mat(DLV(i),CD_ids(j),t) = RP_mat(DLV(i),CD_ids(j)).*mean(pastLB_ratestore(CD_index(j,:)));
                         
                         adj_lender_lend_mat(DLV(i),CD_ids(j),t) = (IBrate_mat(DLV(i),CD_ids(j),t))...
                            .*lender_lend_mat(DLV(i),CD_ids(j));
                        
                        pastLB_Mrel_t      = sprintf(' %.3g ',pastLB_tstore(CD_index(j,:)));
                        pastLB_Mrel_loan   = sprintf(' %.3g ',pastLB_loanstore(CD_index(j,:)));
                        pastLB_Mrel_rep    = sprintf(' %.3g ',pastLB_repstore(CD_index(j,:)));
                        pastLB_Mrel_ratio  = sprintf(' %.3g ',(pastLB_repstore(CD_index(j,:))./pastLB_loanstore(CD_index(j,:)))*100);
                        
                        myprint(writeoption,fileID_D,'Lender %d and borrower %d had MULTIPLE past L/B relationships in periods: [%s]\r\n',...
                        DLV(i),CD_ids(j),pastLB_Mrel_t);
                    
                            myprint(writeoption,fileID_D,'--> Past loans: %s\r\n',pastLB_Mrel_loan);
                            myprint(writeoption,fileID_D,'--> Past repayment: %s\r\n',pastLB_Mrel_rep);
                            myprint(writeoption,fileID_D,'--> Loan to repayment ratio: %s\r\n',pastLB_Mrel_ratio);
                            myprint(writeoption,fileID_D,'--> Sets a counterparty risk premium of %.3f in period %d\r\n',...
                                RP_mat(DLV(i),CD_ids(j))-1,t);
                            myprint(writeoption,fileID_D,'==> The interbank rate for this transaction is %.2f percent\r\n',...
                                (IBrate_mat(DLV(i),CD_ids(j),t)-1)*100);
                            
                            myprint(writeoption,fileID_D,'-------------------------------------------------------------------------------------------\r\n');
                        
                        CD_index = [];
                        %noCD_index = [];
                        pastLB_Mrel_t = [];
                                  
                    end
 
                    for j = 1:length(noCD_ids) % Borrowers with unique past L-B relationship with current lender
                        
                        noCD_index(j,:) = find(pastLB_Bstore == noCD_ids(j));
                        
                        RP_mat(DLV(i),noCD_ids(j)) = ones(1,numel(pastLB_tstore(noCD_index(j,:))))...
                            +(pastLB_loanstore(noCD_index(j,:))-pastLB_repstore(noCD_index(j,:)))...
                            .*(pastLB_loanstoreav(noCD_index(j,:))/pastLB_repstoreav(noCD_index(j,:)));
                        
                        IBrate_mat(DLV(i),noCD_ids(j),t) = RP_mat(DLV(i),noCD_ids(j)).*pastLB_ratestore(noCD_index(j,:));
                        
                        adj_lender_lend_mat(DLV(i),noCD_ids(j),t) = (IBrate_mat(DLV(i),noCD_ids(j),t))...
                            .*lender_lend_mat(DLV(i),noCD_ids(j));
                                                
                        myprint(writeoption,fileID_D,'Lender %d and borrower %d had ONE past L/B relationship in period %d\r\n',...
                            DLV(i),noCD_ids(j),pastLB_tstore(noCD_index(j,:)));
                        
                            myprint(writeoption,fileID_D,'--> Expected repayment: %.3f\r\n',pastLB_loanstore(noCD_index(j,:)));
                            myprint(writeoption,fileID_D,'--> Actual repayment: %.3f\r\n',pastLB_repstore(noCD_index(j,:)));
                            myprint(writeoption,fileID_D,'--> Repayment rate: %d percent\r\n',(pastLB_repstore(noCD_index(j,:))./pastLB_loanstore(noCD_index(j,:)))*100);
                            myprint(writeoption,fileID_D,'--> Sets a counterparty risk premium of %.2f in period %d\r\n',...
                                RP_mat(DLV(i),noCD_ids(j))-1,t);
                            myprint(writeoption,fileID_D,'==> The interbank rate for this transaction is %.2f percent\r\n',...
                                (IBrate_mat(DLV(i),noCD_ids(j),t)-1)*100);
                    
                     myprint(writeoption,fileID_D,'-------------------------------------------------------------------------------------------\r\n');
                        
                        %CD_index   = [];
                        noCD_index = [];

                    end
                    
                else % No duplicate values (multiple past L-B relationships) present amongst current borrowers
                    
                    RP_mat(DLV(i),pastLB_Bstore) = ones(1,numel(pastLB_Bstore))...
                        +(pastLB_loanstore-pastLB_repstore).*(pastLB_loanstoreav/pastLB_repstoreav);
                
                    IBrate_mat(DLV(i),pastLB_Bstore,t) = RP_mat(DLV(i),pastLB_Bstore).*pastLB_ratestore;
                    
                    adj_lender_lend_mat(DLV(i),pastLB_Bstore,t) = IBrate_mat(DLV(i),pastLB_Bstore,t)...
                        .*lender_lend_mat(DLV(i),pastLB_Bstore);
                    
                    %CD_index   = [];
                    %noCD_index = [];
                
                end
            end
                    
            % Identify borrowers for which the STOP condition is not satisfied (no past L-B relationship)
                
                nopastLB_Bstore = setdiff(banks(DLV(i)).borrowing_cps,pastLB_Bstore);
            
                RP_mat(DLV(i),nopastLB_Bstore)       = ones(1,numel(nopastLB_Bstore));
                IBrate_mat(DLV(i),nopastLB_Bstore,t) = (1+r_b).*ones(1,numel(nopastLB_Bstore));
                
                adj_lender_lend_mat(DLV(i),nopastLB_Bstore,t) = IBrate_mat(DLV(i),nopastLB_Bstore,t)...
                    .*lender_lend_mat(DLV(i),nopastLB_Bstore);
                        
                for j=1:length(nopastLB_Bstore)
                    
                     myprint(writeoption,fileID_D,'Lender %d and borrower %d did NOT have a past L/B relationship\r\n',...
                        DLV(i),nopastLB_Bstore(j));
                
                        myprint(writeoption,fileID_D,'--> Sets the default counterparty risk premium of %.2f  in period %d\r\n',...
                            RP_mat(DLV(i),nopastLB_Bstore(j))-1,t);
                
                        myprint(writeoption,fileID_D,'==> The interbank rate for this transaction is %.2f percent\r\n',...
                            (IBrate_mat(DLV(i),nopastLB_Bstore(j),t)-1)*100);
                                    
                end                  
        end  
        
        banks(DLV(i)).IBM.L_bil_exp_repay    = adj_lender_lend_mat(DLV(i),banks(DLV(i)).borrowing_cps,t);
        banks(DLV(i)).IBM.L_tot_exp_repay(t) = sum(adj_lender_lend_mat(DLV(i),banks(DLV(i)).borrowing_cps,t));
        
    end  
       myprint(writeoption,fileID_D,'-------------------------------------------------------------------------------------------\r\n');      
end

%-------------------------------------------------------------------------
%% BALANCE SHEET UPDATING
%-------------------------------------------------------------------------

%current_IBparticipants    = sort([DLV' DBV']);
%current_IBnonparticipants = setdiff(ActiveBanks,current_IBparticipants);

for i = ActiveBanks
        
    if ismember(i,DLV)
        banks(i).balancesheet.liabilities.IB_borrowing(t)   = 0;
        banks(i).balancesheet.assets.IB_lending(t)          = banks(i).IBM.L_tot_loans(t);
    elseif ismember(i,DBV)
        banks(i).balancesheet.assets.IB_lending(t)          = 0;
        banks(i).balancesheet.liabilities.IB_borrowing(t)   =  banks(i).IBM.B_tot_loans(t);
   % elseif ismember(i,B_noIB)
        %banks(i).balancesheet.assets.IB_lending(t)          = 0;
        %banks(i).balancesheet.liabilities.IB_borrowing(t)   = 0; 
    %elseif ismember(i,L_noIB)
      % banks(i).balancesheet.assets.IB_lending(t)          = 0;
         %banks(i).balancesheet.liabilities.IB_borrowing(t)   = 0;   
    end

% Assets

banks(i).balancesheet.assets.external_assets(t,tau) = banks(i).balancesheet.assets.external_assets(t,tau-1);
   % fprintf('error with bank %d\n',i);
    banks(i).balancesheet.assets.total(t,tau) = ...
          banks(i).balancesheet.assets.cash(t,tau)...
        + banks(i).balancesheet.assets.investment(t)...
        + banks(i).balancesheet.assets.external_assets(t,tau)...
        + banks(i).balancesheet.assets.IB_lending(t); 
    
% Liabilities    
    
    banks(i).balancesheet.liabilities.capital(t,tau) =  banks(i).balancesheet.assets.total(t,tau)...
        -(banks(i).balancesheet.liabilities.deposits(t,tau)...
        +banks(i).balancesheet.liabilities.IB_borrowing(t));
    
    banks(i).balancesheet.liabilities.total(t,tau) = banks(i).balancesheet.assets.total(t,tau);
end

%-------------------------------------------------------------------------
%% COLLECTING OUTPUT
%-------------------------------------------------------------------------

NT_matrices(:,t,1) =  shocksign_mat(:,t);
NT_matrices(:,t,2) =  designated_lender_mat';
NT_matrices(:,t,3) =  designated_borr_mat';

NNT_matrices(:,:,t,1) =  expost_adjmat;
NNT_matrices(:,:,t,2) = borr_req_mat(:,:,t);   
NNT_matrices(:,:,t,3) = borr_rec_mat;  
NNT_matrices(:,:,t,4) = borr_rep_mat(:,:,t); 
NNT_matrices(:,:,t,5) = lender_lend_mat;   
NNT_matrices(:,:,t,6) = adj_lender_lend_mat(:,:,t);
NNT_matrices(:,:,t,7) = RP_mat;
NNT_matrices(:,:,t,8) = IBrate_mat(:,:,t);  

%clearvars shocksign_mat designated_lender_mat designated_borr_mat...
   % expost_adjmat borr_req_mat borr_rec_mat borr_rep_mat lender_lend_mat...
    %adj_lender_lend_mat RP_mat IBrate_mat            
                       
end
