function Results_IBM(fig_output,simtype,Results,CB_TOTallotment,T_sim,shocktime,FRFAtime,...
    results_IBM_format,includebounds,outputlabel)

if strcmp(simtype,'baseline')
    fig_output_IBM = strcat(fig_output,'Results/baseline/IBM/');
elseif  strcmp(simtype,'crisis')
    fig_output_IBM = strcat(fig_output,'Results/crisis/IBM/');
end

%Results_agg = Results(:,:,4);

if strcmp(results_IBM_format,'agg_')
    Results_plot =  Results(:,:,4,:);
elseif strcmp(results_IBM_format,'av_')
    Results_plot    = Results(:,:,1,:); % Plot average values
    Results_plot_LB = Results(:,:,2,:); % Plot max and min around average
    Results_plot_UB = Results(:,:,3,:);
end
%-----------------------------------------------------
% Phase 1 Volumes
%-----------------------------------------------------

% Plotting requests and loans separately

figure
subplot(1,2,1) % Interbank loan requests vs actual loans (Baseline)
    if strcmp(results_IBM_format,'agg_')
        plot(Results_plot(8,:),'Color','r','LineWidth',1.1) % Total Loan requests
        hold on
        plot(Results_plot(11,:),'Color','b','LineWidth',1.1) % Total Loans
        title('Interbank loan volumes')
        hold off
    elseif strcmp(results_IBM_format,'av_')
        if strcmp(includebounds,'BOn')
            boundedline(1:T_sim,Results_plot(8,:),[Results_plot_LB(8,:);Results_plot_UB(8,:)]','alpha')
            hold on
            boundedline(1:T_sim,Results_plot(11,:),[Results_plot_LB(11,:);Results_plot_UB(10,:)]','alpha');
        elseif strcmp(includebounds,'BOff')
            if strcmp(simtype,'crisis')
                rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),2],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
                hold on;
                %plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
                %hold on;
            end
            h1 = plot(Results_plot(8,:,1),'Color','r','LineWidth',1.1);
            hold on
            h2 = plot(Results_plot(11,:,1),'Color','b','LineWidth',1.1);
            xlim([0 T_sim]);
        end
    end
    grid on;
    if strcmp(simtype,'baseline')
        title({'(a) Loan volumes','CB refinancing OFF'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(a) FRFA OFF'},'Interpreter','latex')
    end
    legend([h1 h2],{'Requests','Loans'},'Location','best','FontSize',6,'Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
subplot(1,2,2) % Interbank loan requests vs actual loans (with FRFA)
    if strcmp(results_IBM_format,'agg_')
        plot(Results_plot(8,:,2),'Color','b','LineWidth',1.1) % Total Loan requests
        hold on
        plot(Results_plot(11,:,2),'Color','r','LineWidth',1.1) % Total Loans
        title('Interbank loan volumes')
        hold off
    elseif strcmp(results_IBM_format,'av_')
        if strcmp(includebounds,'BOn')
            boundedline(1:T_sim,Results_plot(8,:),[Results_plot_LB(8,:);Results_plot_UB(8,:)]','alpha')
            hold on
            boundedline(1:T_sim,Results_plot(11,:),[Results_plot_LB(11,:);Results_plot_UB(11,:)]','alpha');
        elseif strcmp(includebounds,'BOff')
            if strcmp(simtype,'crisis')
                rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),2],...
                    'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
                hold on;
                plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
                hold on;
            end
            h1 = plot(Results_plot(8,:,2),'Color','r','LineWidth',1.1);
            hold on
            h2 = plot(Results_plot(11,:,2),'Color','b','LineWidth',1.1);
            xlim([0 T_sim]);
        end
    end
    grid on;
    if strcmp(simtype,'baseline')
        title({'(a) Loan volumes','CB refinancing ON'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(b) FRFA ON'},'Interpreter','latex')
    end
    %legend([h1 h2],{'Requests','Loans'},'Location','best','FontSize',8,'Interpreter','latex') 
grid on;
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBM_phase1_vol_',outputlabel,'.pdf'));
    
% Plotting difference between requests and loans in one line     
figure
if strcmp(results_IBM_format,'agg_')
    plot(Results_plot(9,:),'Color','b','LineWidth',1.1) % Total Loan requests
    title('Lender hoarding')
elseif strcmp(results_IBM_format,'av_')
    if strcmp(includebounds,'BOn')
        boundedline(1:T_sim,Results_plot(9,:),[Results_plot_LB(9,:);Results_plot_UB(9,:)]','alpha');
    elseif strcmp(includebounds,'BOff')
        if strcmp(simtype,'crisis')
            rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),1.8],...
                'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
            hold on;
            plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
            hold on;
        end
        h1 = plot(Results_plot(9,:,1),'Color','r','LineWidth',1.1);
        hold on
        h2 = plot(Results_plot(9,:,2),'Color','b','LineWidth',1.1);
        xlim([0 T_sim]);
    end
end
title('Lender hoarding','Interpreter','latex')
if strcmp(simtype,'baseline')
    legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},...
        'Location','best','FontSize',6,'Interpreter','latex')
elseif strcmp(simtype,'crisis')
    legend([h1 h2],{'FRFA OFF','FRFA ON'},...
        'Location','best','FontSize',6,'Interpreter','latex')
end
grid on;
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBM_phase1_diff_',outputlabel,'.pdf'));


% Decomposing lender hoarding
figure
subplot(1,2,1)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),0.6],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on;
    end
    bar([(Results_plot(8,:,1)-Results_plot(9,:,1))',(Results_plot(9,:,1)-Results_plot(11,:,1))'],'stacked')
    xlim([0 T_sim]);
    if strcmp(simtype,'baseline')
        title({'(a) Loan decomposition','CB refinancing ON'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(a) FRFA OFF'},'Interpreter','latex')
    end
    grid on;
    legend({'Liquidity motive','Precautionary motive'},'Location','best','FontSize',6,'Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
subplot(1,2,2)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),0.5],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on;
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    bar([(Results_plot(8,:,2)-Results_plot(9,:,2))',(Results_plot(9,:,2)-Results_plot(11,:,2))'],'stacked')
    xlim([0 T_sim]);
    if strcmp(simtype,'baseline')
        title({'(a) Loan decomposition','CB refinancing ON'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(b) FRFA ON'},'Interpreter','latex')
    end
    grid on;
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 18, 8],'PaperUnits','Centimeters','PaperSize',[18,8])
print(gcf,'-dpdf',strcat(fig_output_IBM,'loan_decomp_',outputlabel,'.pdf'));


%-----------------------------------------------------
% Phase 2 Volumes
%-----------------------------------------------------

% Plotting expected and final loan repayment separately
figure
subplot(1,2,1) % Expected vs. actual interbank loan repayment (baseline)
    if strcmp(results_IBM_format,'agg_')
        plot(Results_plot(11,:),'Color','b','LineWidth',1.1) % Total expected repayment
        hold on
        plot(Results_plot(12,:),'Color','k','LineWidth',1.1) % Total Loans
        title('Interbank loan volumes')
        hold off
    elseif strcmp(results_IBM_format,'av_')
        if strcmp(includebounds,'BOn')
            boundedline(1:T_sim,Results_plot(8,:),[Results_plot_LB(8,:);Results_plot_UB(8,:)]','alpha')
            hold on
            boundedline(1:T_sim,Results_plot(10,:),[Results_plot_LB(10,:);Results_plot_UB(10,:)]','alpha');
        elseif strcmp(includebounds,'BOff')
            if strcmp(simtype,'crisis')
                rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),1.5],...
                    'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
                hold on;
                %plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
                %hold on;
            end
            h1 = plot(Results_plot(12,:,1),'Color','r','LineWidth',1.1);
            hold on
            h2 = plot(Results_plot(13,:,1),'Color','b','LineWidth',1.1);
            xlim([0 T_sim]);
        end
    end
    grid on;
    if strcmp(simtype,'baseline')
        title({'(a) Repayment volumes','CB refinancing OFF'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(a) FRFA OFF'},'Interpreter','latex')
    end
    legend([h1 h2],{'Expected repayment','Final repayment'},...
        'Location','best','FontSize',6,'Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
subplot(1,2,2) % Expected vs. actual interbank loan repayment (FRFA active)
    if strcmp(results_IBM_format,'agg_')
        plot(Results_plot(12,:,2),'Color','r','LineWidth',1.1) % Total Loan requests
        hold on
        plot(Results_plot(13,:,2),'Color','b','LineWidth',1.1) % Total Loans
        title('Interbank repayment volumes')
        hold off
    elseif strcmp(results_IBM_format,'av_')
        if strcmp(includebounds,'BOn')
            boundedline(1:T_sim,Results_plot(8,:),[Results_plot_LB(8,:);Results_plot_UB(8,:)]','alpha')
            hold on
            boundedline(1:T_sim,Results_plot(10,:),[Results_plot_LB(10,:);Results_plot_UB(10,:)]','alpha');
        elseif strcmp(includebounds,'BOff')
            if strcmp(simtype,'crisis')
                rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),2],'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
                hold on;
                plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',1.5);
                hold on;
            end
            h1 = plot(Results_plot(12,:,2),'Color','r','LineWidth',1.1);
            hold on
            h2 = plot(Results_plot(13,:,2),'Color','b','LineWidth',1.1);
            xlim([0 T_sim]);
        end
    end
    grid on;
    if strcmp(simtype,'baseline')
        title({'(a) Repayment volumes','CB refinancing ON'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(b) FRFA ON'},'Interpreter','latex')
    end
    %legend([h1 h2],{'Requests','Loans'},'Location','best','FontSize',8,'Interpreter','latex')  
 grid on;
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBM_phase2_vol_',outputlabel,'.pdf'));  
    
% Plotting difference between expected and final repayment in one line
figure
if strcmp(results_IBM_format,'agg_')
    plot(Results_plot(9,:),'Color','b','LineWidth',1.1) % Total Loan requests
    title('Lender hoarding')
elseif strcmp(results_IBM_format,'av_')
    if strcmp(includebounds,'BOn')
        boundedline(1:T_sim,Results_plot(9,:),[Results_plot_LB(9,:);Results_plot_UB(9,:)]','alpha');
    elseif strcmp(includebounds,'BOff')
        if strcmp(simtype,'crisis')
            rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),0.04],...
            'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
             hold on;
             plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
             hold on;
             xlim([0 T_sim]);
        end
        h1 = plot(Results_plot(12,:,1)-Results_plot(13,:,1),'Color','r','LineWidth',1.1);
        hold on
        h2 = plot(Results_plot(12,:,2)-Results_plot(13,:,2),'Color','b','LineWidth',1.1);
     end
end
%title('(c) Borrower defaults','Interpreter','latex')
if strcmp(simtype,'baseline')
    legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},'Location','best','FontSize',8,'Interpreter','latex')
elseif strcmp(simtype,'crisis')
    legend([h1 h2],{'FRFA OFF','FRFA ON'},'Location','best','FontSize',8,'Interpreter','latex')
end
grid on;
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBM_phase2_diff_',outputlabel,'.pdf'));


%-----------------------------------------------------
% Interbank Rate
%-----------------------------------------------------

figure
subplot(1,2,1)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),4,shocktime(end)-shocktime(1),2],...
            'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
    end
    h1 = plot(1:T_sim,(Results_plot(16,:,1)-1)*100);
    hold on;
    h2 = plot(1:T_sim,(Results_plot_UB(16,:,1)-1)*100);
    if strcmp(simtype,'baseline')
        title({'(a) Interbank rate','CB refinancing OFF'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(a) FRFA OFF'},'Interpreter','latex')
    end
    xlim([0 T_sim]);
    grid on;
    legend([h1 h2],{'$\bar{r}^{b}$','$\max\left({r}^{b}\right)$'},'Location','best','FontSize',8,'Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
subplot(1,2,2)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),4,shocktime(end)-shocktime(1),2],...
            'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(1:T_sim,(Results_plot(16,:,2)-1)*100);
    hold on;
    h2 = plot(1:T_sim,(Results_plot_UB(16,:,2)-1)*100);
    if strcmp(simtype,'baseline')
        title({'(a) Interbank rate','CB refinancing ON'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(b) FRFA ON'},'Interpreter','latex')
    end
    xlim([0 T_sim]);  
    grid on;
%----------------------------------------------------------------------------------------------------------
set(gcf,'renderer','painters');
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 18, 8],'PaperUnits','Centimeters','PaperSize',[18,8])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBrates_',outputlabel,'.pdf'));

%--------------------------------------------------
% Central bank liquidity provision
%--------------------------------------------------

figure
if strcmp(includebounds,'BOn')
    confplot(1:T_sim,CB_TOTallotment(1,:),CB_TOTallotment(2,:),CB_TOTallotment(3,:))
elseif strcmp(includebounds,'BOff')
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),60],...
            'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on;
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(CB_TOTallotment(1,:,2),'LineWidth',1.1);
    hold on
    h2 = plot(CB_TOTallotment(2,:,2),'LineWidth',1.1);
    if strcmp(simtype,'crisis')
        h3 =plot(CB_TOTallotment(3,:,2),'LineWidth',1.1);
    end
    xlim([0 T_sim]);
end
grid on;
if strcmp(simtype,'baseline')
    legend([h1 h2],{'Borrower refinancing','Lender refinancing'},...
        'Location','best','FontSize',6,'Interpreter','latex');
elseif strcmp(simtype,'crisis')
     legend([h1 h2 h3],{'Borrower refinancing','Lender refinancing','FRFA refinancing'},...
        'Location','best','FontSize',6,'Interpreter','latex');
end
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'CBliquidity_',outputlabel,'.pdf'));



end