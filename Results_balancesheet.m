function Results_balancesheet(fig_output,Results,T,T_sim,normalizegraph,results_BS_format,includebounds,outputlabel,allsimlabels)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure output: Balance sheet variables evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results_av  = Results(:,:,1);
Results_min = Results(:,:,2);
Results_max = Results(:,:,3);
Results_agg = Results(:,:,4);

fig_output_BS = strcat(fig_output,'Results/balancesheet/');

 % Select whether to normalize results by dividing by total assets

if strcmp(normalizegraph,'NormOn_')
    normfactor_agg = Results_agg(1,:);
    normfactor_av  = Results_av(1,:);
    normfactor_min = Results_min(1,:);
    normfactor_max = Results_max(1,:);
elseif strcmp(normalizegraph,'NormOff_')
    normfactor_agg = ones(1,T);
    normfactor_av  = ones(1,T);
    normfactor_min = ones(1,T);
    normfactor_max = ones(1,T);
end

if strcmp(results_BS_format,'agg_')
    Results_plot = Results_agg;
elseif strcmp(results_BS_format,'av_')
    Results_plot    = Results_av; % Plot average values
    Results_plot_LB = Results_min; % Plot max and min around average
    Results_plot_UB = Results_max;
end

%--------------------------------------------------
% Aggregate balance sheet size
%--------------------------------------------------

figure
if strcmp(results_BS_format,'agg_')
    plot(Results_plot(1,:),'LineWidth',1.1)
elseif strcmp(results_BS_format,'av_')
    if strcmp(includebounds,'BOn')
        confplot(1:T_sim,Results_plot(1,:),Results_plot_LB(1,:),Results_plot_UB(1,:))
    elseif strcmp(includebounds,'BOff')
        plot(Results_plot(1,:))
    end
end
grid on;
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_BS,'totalassets_',outputlabel,'_',allsimlabels,'.pdf'));

%--------------------------------------------------
% Evolution of balance sheet componenents
%--------------------------------------------------

title_vec = {'Cash','External assets','Deposits','Capital'};
splot_vec = [2 3 6 7];

figure
for k = 1:4
subplot(2,2,k)
    if strcmp(results_BS_format,'agg_')
        plot(Results_plot(splot_vec(k),:)./normfactor_agg,'LineWidth',1.1)
    elseif strcmp(results_BS_format,'av_')
        if strcmp(includebounds,'BOn')
        confplot(1:T_sim,Results_plot(splot_vec(k),:)./normfactor_av,...
                 Results_plot_LB(splot_vec(k),:)./normfactor_min,...
                 Results_plot_UB(splot_vec(k),:)./normfactor_max)
        elseif  strcmp(includebounds,'BOff')
            plot(Results_plot(splot_vec(k),:)./normfactor_av)
        end
    end
    grid on;
    title(title_vec(k),'Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
end
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_BS,'balancesheet_',outputlabel,'_',allsimlabels,'.pdf'));

%--------------------------------------------------
% Investment
%--------------------------------------------------

figure
plot(Results_plot(4,:)./normfactor_av,'Color','b','LineWidth',1.1)
hold on
plot(Results_plot(5,:)./normfactor_av,'Color','r','LineWidth',1.1)
grid on;
%grid minor;
legend({'Total planned investment','Final investment'},'Location','best','FontSize',8,'Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_BS,'investment_',outputlabel,'_',allsimlabels,'.pdf'));


end
