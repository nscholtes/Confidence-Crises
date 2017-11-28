function Results_IBM(fig_output,Results,CB_TOTallotment,T_sim,results_IBM_format,includebounds,outputlabel,allsimlabels)

Results_av  = Results(:,:,1);
Results_min = Results(:,:,2);
Results_max = Results(:,:,3);
Results_agg = Results(:,:,4);

fig_output_IBM = strcat(fig_output,'Results/IBM/');

if strcmp(results_IBM_format,'agg_')
    Results_plot = Results_agg;
elseif strcmp(results_IBM_format,'av_')
    Results_plot    = Results_av; % Plot average values
    Results_plot_LB = Results_min; % Plot max and min around average
    Results_plot_UB = Results_max;
end

%-----------------------------------------------------
% Phase 1 Volumes
%-----------------------------------------------------

figure
subplot(1,2,1) % Interbank loan requests vs actual loans
    if strcmp(results_IBM_format,'agg_')
        plot(Results_plot(8,:),'Color','b','LineWidth',1.1) % Total Loan requests
        hold on
        plot(Results_plot(10,:),'Color','k','LineWidth',1.1) % Total Loans
        title('Interbank loan volumes')
        hold off
    elseif strcmp(results_IBM_format,'av_')
        if strcmp(includebounds,'BOn')
            boundedline(1:T_sim,Results_plot(8,:),[Results_plot_LB(8,:);Results_plot_UB(8,:)]','alpha')
            hold on
            boundedline(1:T_sim,Results_plot(10,:),[Results_plot_LB(10,:);Results_plot_UB(10,:)]','alpha');
        elseif strcmp(includebounds,'BOff')
            plot(Results_plot(8,:),'LineWidth',1.1)
            hold on
            plot(Results_plot(10,:),'LineWidth',1.1)
        end
    end
    grid on;
    title('(a) Interbank loan volumes','Interpreter','latex')
    legend({'Requests','Loans'},'Location','best','FontSize',8,'Interpreter','latex')
subplot(1,2,2) % Hoarding by lenders
    if strcmp(results_IBM_format,'agg_')
        plot(Results_plot(9,:),'Color','b','LineWidth',1.1) % Total Loan requests
        title('Lender hoarding')
    elseif strcmp(results_IBM_format,'av_')
       if strcmp(includebounds,'BOn')
            boundedline(1:T_sim,Results_plot(9,:),[Results_plot_LB(9,:);Results_plot_UB(9,:)]','alpha');
       elseif strcmp(includebounds,'BOff')
           plot(Results_plot(9,:))
       end
    end
title('(b) Lender hoarding','Interpreter','latex')
grid on;
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBM_phase1_',outputlabel,'_',allsimlabels,'.pdf'));

%-----------------------------------------------------
% Phase 2 Volumes
%-----------------------------------------------------

figure
subplot(1,2,1) % Expected interbank loan repayment vs. Actual repayment
    if strcmp(results_IBM_format,'agg_')
        plot(Results_plot(11,:),'Color','b','LineWidth',1.1) % Total Loan requests
        hold on
        plot(Results_plot(12,:),'Color','r','LineWidth',1.1) % Total Loans
        hold off
    elseif strcmp(results_IBM_format,'av_')
        if strcmp(includebounds,'BOn')
            boundedline(1:T_sim,Results_plot(11,:),[Results_plot_LB(11,:);Results_plot_UB(11,:)]','alpha',...
                        1:T_sim,Results_plot(12,:),[Results_plot_LB(12,:);Results_plot_UB(12,:)]','alpha');
        elseif strcmp(includebounds,'BOff')
            plot(Results_plot(11,:),'b','LineWidth',1.1)
            hold on
            plot(Results_plot(12,:),'r','LineWidth',1.1)
        end
    end
    title('(a) Interbank repayment volumes','Interpreter','latex')
    grid on;
    legend({'Requests','Loans'},'Location','best','FontSize',8,'Interpreter','latex')
subplot(1,2,2) % Hoarding by lenders
    if strcmp(results_IBM_format,'agg_')
        plot(Results_plot(11,:)-Results_plot(12,:),'Color','b','LineWidth',1.1) % Total Loan requests
    elseif strcmp(results_IBM_format,'av_')
        if strcmp(includebounds,'BOn')
        boundedline(1:T_sim,Results_plot(11,:)-Results_plot(12,:),...
            [Results_plot_LB(11,:)-Results_plot_LB(12,:);Results_plot_UB(11,:)-Results_plot_UB(12,:)]','alpha');
        elseif strcmp(includebounds,'BOff')
            plot((Results_plot(11,:)-Results_plot(12,:)),'b','LineWidth',1.1)
        end
    end
title('(b) Loan defaults','Interpreter','latex')
grid on; 
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBM_phase2_',outputlabel,'_',allsimlabels,'.pdf'));



%-----------------------------------------------------
% Interbank Rate
%-----------------------------------------------------

figure
if strcmp(includebounds,'BOn')
    confplot(1:T_sim,Results_av(15,:),Results_min(15,:),Results_max(15,:))
elseif strcmp(includebounds,'BOff')
    plot(Results_av(15,:),'Color','k','LineWidth',1.1)
end
xlabel('Iteration step','Interpreter','latex')
grid on;
%legend({'$\bar{r}^{b,t}$''},'Location','best','FontSize',8,'Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBrates_',outputlabel,'_',allsimlabels,'.pdf'));

%--------------------------------------------------
% Central bank liquidity provision
%--------------------------------------------------

figure
if strcmp(includebounds,'BOn')
    confplot(1:T_sim,CB_TOTallotment(1,:),CB_TOTallotment(2,:),CB_TOTallotment(3,:))
elseif strcmp(includebounds,'BOff')
    plot(CB_TOTallotment(1,:))
    hold on
    plot(CB_TOTallotment(2,:))
end
grid on;
legend({'Borrower refinancing','Lender refinancing'},'Location','best','FontSize',8,'Interpreter','latex');
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'CBliquidity_',outputlabel,'_',allsimlabels,'.pdf'));


end