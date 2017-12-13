function Results_securitiesmarket(fig_output,simtype,Results,assetprices,T,T_sim,shocktime,FRFAtime,...
    results_SM_format,includebounds,outputlabel)

if strcmp(simtype,'baseline')
    fig_output_SM = strcat(fig_output,'Results/baseline/securitiesmarket/');
elseif  strcmp(simtype,'crisis')
    fig_output_SM = strcat(fig_output,'Results/crisis/securitiesmarket/');
end

if strcmp(results_SM_format,'agg_')
    Results_plot =  Results(:,:,4,:);
elseif strcmp(results_SM_format,'av_')
    Results_plot    = Results(:,:,1,:); % Plot average values
    Results_plot_LB = Results(:,:,2,:); % Plot max and min around average
    Results_plot_UB = Results(:,:,3,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Asset prices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AP_axis = cell(1,2*T);

deltaAP_APS = zeros(1,T);
deltaAP_MIF = zeros(1,T);

%for t=1:T    
 %   AP_axis{(2*t)-1} =  strcat('APS_','{',num2str(t),'}');
 %   AP_axis{2*t}     =  strcat('MIF_','{',num2str(t),'}');
%end    

%AP_axis = ['0',AP_axis];

%x_axis = [1 linspace(2,4*T,2*T)];
%keep_assetprices = assetprices(x_axis,:);

APSindices = linspace(2,(4*T)-2,T);
MIFindices = linspace(4,4*T,T);

for t = 1:T
    for S = 1:2
    if t>1
        deltaAP_APS(S,t) = nanmean(((assetprices(APSindices(t),:,S) - assetprices(MIFindices(t-1),:,S))/(assetprices(MIFindices(t-1),:,S)))*100);
    end
    deltaAP_MIF(S,t) = nanmean(((assetprices(MIFindices(t),:,S) - assetprices(APSindices(t),:,S))/(assetprices(APSindices(t),:,S)))*100);

    end
end

% Average asset prices over simulation period
figure
if strcmp(simtype,'crisis')
    rectangle('Position',[4*shocktime(1),0,4*(shocktime(end)-shocktime(1)),1.2],...
    'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
    hold on
    plot([4*FRFAtime(1), 4*FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
    hold on;
end
    h1 = plot(mean(assetprices(:,:,1),2),'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(mean(assetprices(:,:,2),2),'Color','b','LineWidth',1.1);
    xlabel('Iteration step','Interpreter','latex')
    grid on;
    xlim([0 4*T_sim]);
    if strcmp(simtype,'baseline')
        legend([h1 h2],{'$\bar{p}$: CB refinancing OFF','$\bar{p}$: CB refinancing ON'},...
            'Location','best','FontSize',6,'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        legend([h1 h2],{'$\bar{p}$: FRFA OFF','$\bar{p}$: FRFA ON'},...
            'Location','best','FontSize',6,'Interpreter','latex')
    end
    %set(gca,'XLim',[x_axis(1) x_axis(end)+1],'XTick',x_axis,'XTickLabel',AP_axis,'XTickLabelRotation',45)
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_SM,'assetprices','.pdf'));

% Asset prices changes due to exogenous shocks and endogenous firesales
figure
subplot(1,2,1)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),-10,shocktime(end)-shocktime(1),12],...
            'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        %plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        %hold on;
    end
    h1 = plot(deltaAP_APS(1,:),'Color','b','LineWidth',1.1);
    hold on
    h2 = plot(deltaAP_MIF(1,:),'Color','r','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on
    if strcmp(simtype,'baseline')
        title({'(a) Asset price dynamics','CB refinancing OFF'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(b) Asset price dynamics','FRFA OFF'},'Interpreter','latex')
    end
    legend([h1 h2],{'$\Delta\bar{p}^{APS}$','$\Delta\bar{p}^{MIF}$'},...
        'Location','best','FontSize',6,'Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
subplot(1,2,2)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),-10,shocktime(end)-shocktime(1),12],...
            'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(deltaAP_APS(2,:),'Color','b','LineWidth',1.1);
    hold on
    h2 = plot(deltaAP_MIF(2,:),'Color','r','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on
    if strcmp(simtype,'baseline')
        title({'(a) Asset price dynamics','CB refinancing ON'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(b) Asset price dynamics','FRFA ON'},'Interpreter','latex')
    end
    %legend([h1 h2],{'$\Delta\bar{p}^{APS}$','$\Delta\bar{p}^{MIF}$'},'Location','best','FontSize',8,'Interpreter','latex')
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_SM,'assetprices_components','.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Firesales %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
%----------------------------------------------------------------------------------------------------------
subplot(1,2,1) % Desired and expected firesales: CB refinancing OFF/FRFA OFF
    if strcmp(results_SM_format,'agg_')
        plot(Results_plot(13,:),'Color','b','LineWidth',1.1) % Total expected repayment
        hold on
        plot(Results_plot(14,:),'Color','k','LineWidth',1.1) % Total Loans
        title('Interbank loan volumes')
        hold off
    elseif strcmp(results_SM_format,'av_')
        if strcmp(includebounds,'BOn')
            boundedline(1:T_sim,Results_plot(8,:),[Results_plot_LB(8,:);Results_plot_UB(8,:)]','alpha')
            hold on
            boundedline(1:T_sim,Results_plot(10,:),[Results_plot_LB(10,:);Results_plot_UB(10,:)]','alpha');
        elseif strcmp(includebounds,'BOff')
            if strcmp(simtype,'crisis')
                rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),1.2],...
                    'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
                hold on;
                %plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
                %hold on;
            end
            h1 = plot(Results_plot(13,:,1),'Color','r','LineWidth',1.1);
            hold on
            h2 = plot(Results_plot(14,:,1),'Color','b','LineWidth',1.1); 
            xlim([0 T_sim]);
        end
    end
    grid on;
    if strcmp(simtype,'baseline')
        title({'(a) Firesales','CB refinancing OFF'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(b) Firesales','FRFA OFF'},'Interpreter','latex')
    end
    legend([h1 h2],{'Desired firesales','Final firesales'},'Location','best','FontSize',8,'Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
subplot(1,2,2) % Desired and expected firesales: CB refinancing ON/FRFA ON
    if strcmp(results_SM_format,'agg_')
        plot(Results_plot(13,:,2),'Color','r','LineWidth',1.1) % Total Loan requests
        hold on
        plot(Results_plot(14,:,2),'Color','b','LineWidth',1.1) % Total Loans
        title('Interbank repayment volumes')
        hold off
    elseif strcmp(results_SM_format,'av_')
        if strcmp(includebounds,'BOn')
            boundedline(1:T_sim,Results_plot(8,:),[Results_plot_LB(8,:);Results_plot_UB(8,:)]','alpha')
            hold on
            boundedline(1:T_sim,Results_plot(10,:),[Results_plot_LB(10,:);Results_plot_UB(10,:)]','alpha');
        elseif strcmp(includebounds,'BOff')
            if strcmp(simtype,'crisis')
                rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),1.2],...
                    'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
                hold on;
                plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
                hold on;    
            end
            h1 = plot(Results_plot(13,:,2),'Color','r','LineWidth',1.1);
            hold on
            h2 = plot(Results_plot(14,:,2),'Color','b','LineWidth',1.1);
            xlim([0 T_sim]);
        end
    end
    grid on;
    if strcmp(simtype,'baseline')
        title({'(a) Firesale volumes','CB refinancing ON'},'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title({'(b) Firesale volumes','FRFA ON'},'Interpreter','latex')
    end
    %legend([h1 h2],{'Desired firesales','Final firesales'},'Location','best','FontSize',8,'Interpreter','latex')
xlabel('Iteration step','Interpreter','latex')
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_SM,'firesales_vol_',outputlabel,'.pdf'));

%----------------------------------------------------------------------------------------------------------
% Difference between expected and final firesales
figure
if strcmp(results_SM_format,'agg_')
    plot(Results_plot(9,:),'Color','b','LineWidth',1.1) % Total Loan requests
    title('Lender hoarding')
elseif strcmp(results_SM_format,'av_')
    if strcmp(includebounds,'BOn')
        boundedline(1:T_sim,Results_plot(9,:),[Results_plot_LB(9,:);Results_plot_UB(9,:)]','alpha');
    elseif strcmp(includebounds,'BOff')
        if strcmp(simtype,'crisis')
            rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),0.2],...
                'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
            hold on;
            plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
            hold on;          
        end
        h1 = plot(Results_plot(13,:,1)-Results_plot(14,:,1),'Color','r','LineWidth',1.1);
        hold on
        h2 = plot(Results_plot(13,:,2)-Results_plot(14,:,2),'Color','b','LineWidth',1.1);
        xlim([0 T_sim]);
    end
end
%title('Firesale gap','Interpreter','latex')
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
print(gcf,'-dpdf',strcat(fig_output_SM,'firesales_diff_',outputlabel,'.pdf'));

end