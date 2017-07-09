function Results_securitiesmarket(fig_output,TOT_T_matrices,TM_matrices,total_firesales_vec,assetprices,T)

fig_output_SM = strcat(fig_output,'Results/');

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
    if t>1
        deltaAP_APS(t) = nanmean(((assetprices(APSindices(t),:) - assetprices(MIFindices(t-1),:))/(assetprices(MIFindices(t-1),:)))*100);
    end
    deltaAP_MIF(t) = nanmean(((assetprices(MIFindices(t),:) - assetprices(APSindices(t),:))/(assetprices(APSindices(t),:)))*100);
end

% Average asset prices over simulation period
figure
subplot(1,2,1)
    plot(mean(assetprices,2),'Color','k','LineWidth',1.1)
    xlabel('Iteration step','Interpreter','latex')
    grid on;
    %set(gca,'XLim',[x_axis(1) x_axis(end)+1],'XTick',x_axis,'XTickLabel',AP_axis,'XTickLabelRotation',45)
    %legend({'$\bar{p}^{t}$','$min\{p^{t}\}$','$max\{p^{t}\}$'},'Location','best','FontSize',8,'Interpreter','latex')
subplot(1,2,2)
    plot(deltaAP_APS,'Color','b','LineWidth',1.1)
    hold on
    plot(deltaAP_MIF,'Color','r','LineWidth',1.1)
    grid on
    legend({'$\Delta\bar{p}^{APS}$','$\Delta\bar{p}^{MIF}$'},'Location','best','FontSize',8,'Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_SM,'assetprices.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Firesales %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,2,1)
    plot(TOT_T_matrices(13,:),'Color','b','LineWidth',1.1);
    hold on
    plot(TOT_T_matrices(14,:),'Color','r','LineWidth',1.1);
    grid on; 
    title('(a)','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
    legend({'Desired firesales','Final firesales'},'Location','best','FontSize',8)
subplot(1,2,2)
    plot(TOT_T_matrices(13,:)-TOT_T_matrices(14,:),'Color','k','LineWidth',1.1);
    title('(b)','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_SM,'firesales.pdf'));




end