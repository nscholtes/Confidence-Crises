function Results_balancesheet(fig_output,simtype,Results,T_sim,shocktime,FRFAtime,...
    normalizegraph,results_BS_format,includebounds,outputlabel)

if strcmp(simtype,'baseline')
    fig_output_BS = strcat(fig_output,'Results/baseline/balancesheet/');
elseif  strcmp(simtype,'crisis')
    fig_output_BS = strcat(fig_output,'Results/crisis/balancesheet/');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure output: Balance sheet variables evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Results_av  = Results(:,:,1);
%Results_min = Results(:,:,2);
%Results_max = Results(:,:,3);
%Results_agg = Results(:,:,4);

 % Select whether to normalize results by dividing by total assets

if strcmp(normalizegraph,'NormOn_')
    normfactor_av  = Results(1,:,1,:);
    normfactor_min = Results(1,:,2,:);
    normfactor_max = Results(1,:,3,:);
    normfactor_agg = Results(1,:,4,:);

elseif strcmp(normalizegraph,'NormOff_')
    normfactor_agg = ones(1,T_sim);
    normfactor_av  = ones(1,T_sim);
    normfactor_min = ones(1,T_sim);
    normfactor_max = ones(1,T_sim);
end

if strcmp(results_BS_format,'agg_')
    Results_plot =  Results(:,:,4,:);
elseif strcmp(results_BS_format,'av_')
    Results_plot    = Results(:,:,1,:);  % Plot average values
    Results_plot_LB = Results(:,:,2,:); % Plot max and min around average
    Results_plot_UB = Results(:,:,3,:);
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
        if strcmp(simtype,'crisis')
            rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),100],...
                'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);    
            hold on;
            plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
            hold on;
        end
        h1 = plot(Results_plot(1,:,1),'Color','r','LineWidth',1.1);
        hold on
        h2 = plot(Results_plot(1,:,2),'Color','b','LineWidth',1.1);
        xlim([0 T_sim]);
        if strcmp(simtype,'baseline')
            legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},...
                'Location','best','FontSize',6,'Interpreter','latex')
        elseif strcmp(simtype,'crisis')
            legend([h1 h2],{'FRFA OFF','FRFA ON'},...
                'Location','best','FontSize',6,'Interpreter','latex')
        end
    end
end
grid on;
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_BS,'totalassets_',outputlabel,'.pdf'));

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
            if strcmp(simtype,'crisis')
                rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),50],...
                    'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
                hold on;
                plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
                hold on;
            end
            h1 = plot(Results_plot(splot_vec(k),:,1)./normfactor_av,'Color','r','LineWidth',1.1);
            hold on;
            h2 = plot(Results_plot(splot_vec(k),:,2)./normfactor_av,'Color','b','LineWidth',1.1);
            xlim([0 T_sim]);
            if k == 1
                if strcmp(simtype,'baseline')
                    legend([h1 h2], {'CB refinancing OFF','CB refinancing ON'},...
                        'Location','best','FontSize',6,'Interpreter','latex')
                elseif strcmp(simtype,'crisis')
                    legend([h1 h2],{'FRFA OFF','FRFA ON'},...
                        'Location','best','FontSize',6,'Interpreter','latex')
                end
            end
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
print(gcf,'-dpdf',strcat(fig_output_BS,'balancesheet_',outputlabel,'.pdf'));

%--------------------------------------------------
% Investment
%--------------------------------------------------

figure
subplot(1,2,1)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),1],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on;
        %plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2)
        %hold on;
    end
    h1 = plot(Results_plot(4,:,1)./normfactor_av,'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(Results_plot(5,:,1)./normfactor_av,'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    if strcmp(simtype,'baseline')
        title('(a) CB refinancing OFF','Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title('(a) FRFA OFF','Interpreter','latex')
    end
    legend([h1 h2], {'Total planned investment','Final investment'},...
        'Location','best','FontSize',6,'Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
subplot(1,2,2)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),1],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on;
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2)
        hold on;
    end
    h1 = plot(Results_plot(4,:,2)./normfactor_av,'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(Results_plot(5,:,2)./normfactor_av,'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    if strcmp(simtype,'baseline')
        title('(a) CB refinancing ON','Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title('(a) FRFA ON','Interpreter','latex')
    end
    %legend([h1 h2], {'Total planned investment','Final investment'},'Location','best','FontSize',8,'Interpreter','latex')
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_BS,'investment_',outputlabel,'.pdf'));


end
