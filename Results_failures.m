function [] = Results_failures(fig_output,simtype,n_banks,T_sim,...
    FailCount,cum_Fails,capitalshortfall,N_vec,M_vec,numedges,density,avdegree,mu_A,mu_B,shocktime,FRFAtime)

if strcmp(simtype,'baseline')
    fig_output_F = strcat(fig_output,'Results/baseline/failures/');
elseif  strcmp(simtype,'crisis')
    fig_output_F = strcat(fig_output,'Results/crisis/failures/');
end

% FailCount_av  = FailCount(1,:); FailCount_min = FailCount(2,:); FailCount_max = FailCount(3,:);
% cum_Fails_av  = cum_Fails(1,:); cum_Fails_min = cum_Fails(2,:); cum_Fails_max = cum_Fails(3,:);
% 
% capitalshortfall_av  = capitalshortfall(1,:); capitalshortfal_min = capitalshortfall(2,:); capitalshortfal_max = capitalshortfall(3,:);
% 
% N_av  = N_vec(1,:); N_min = N_vec(2,:); N_max = N_vec(3,:);
% M_av  = M_vec(1,:); M_min = M_vec(2,:); M_max = M_vec(3,:);
% 
% numedges_av = numedges(1,:); numedges_min = numedges(2,:); numedges_max = numedges(3,:);
% density_av  = density(1,:);  density_min  = density(2,:);  density_max  = density(3,:);
% avdegree_av = avdegree(1,:); avdegree_min = avdegree(2,:); avdegree_max = avdegree(3,:);
% 
% mu_A_av = mu_A(1,:); mu_A_min = mu_A(2,:); mu_A_max = mu_A(3,:);
% mu_B_av = mu_B(1,:); mu_B_min = mu_B(2,:); mu_B_max = mu_B(3,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Cumulative number of failures over simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

figure
subplot(1,2,1)
if strcmp(simtype,'crisis')    
    rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),n_banks],...
    'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
    hold on
    plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
    hold on;
end
    h1 = plot(cum_Fails(1,:,1),'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(cum_Fails(1,:,2),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    xlabel('Iteration step','Interpreter','latex')
    grid on;
    title('Cumulative failures','Interpreter','latex')
     if strcmp(simtype,'baseline')
        legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},...
            'Location','northwest','FontSize',8,'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        legend([h1 h2],{'FRFA OFF','FRFA ON'},...
            'Location','northwest','FontSize',8,'Interpreter','latex')
     end
%----------------------------------------------------------------------------------------------------------
subplot(1,2,2)
if strcmp(simtype,'crisis')    
    rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),n_banks],...
    'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
    hold on
    plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
    hold on;
end
    h1 = plot(cumsum(capitalshortfall(1,:,1)),'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(cumsum(capitalshortfall(1,:,2)),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    xlabel('Iteration step','Interpreter','latex')
    grid on;
    title('Cumulative capital shortfall','Interpreter','latex')
     %if strcmp(simtype,'baseline')
        %legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %elseif strcmp(simtype,'crisis')
        %legend([h1 h2],{'FRFA OFF','FRFA ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %end
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_F,'failures+CS','.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Graphing changing network structure due to ABM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%-----------------------------------------------------
% Interbank exposure network
%-----------------------------------------------------

figure
subplot(1,3,1)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),160],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(numedges(1,:,1),'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(numedges(1,:,2),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    title('Number of edges','Interpreter','latex')
    if strcmp(simtype,'baseline')
        legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},...
            'Location','southwest','FontSize',8,'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        legend([h1 h2],{'FRFA OFF','FRFA ON'},...
            'Location','southwest','FontSize',8,'Interpreter','latex')
    end
    xlabel('Iteration step','Interpreter','latex')
%----------------------------------------------------------------------------------------------------------    
subplot(1,3,2)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0.04,shocktime(end)-shocktime(1),0.03],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(density(1,:,1),'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(density(1,:,2),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    title('Network density','Interpreter','latex')
    %if strcmp(simtype,'baseline')
        %legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %elseif strcmp(simtype,'crisis')
        %legend([h1 h2],{'FRFA OFF','FRFA ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %end
    xlabel('Iteration step','Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
subplot(1,3,3)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),5],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(avdegree(1,:,1),'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(avdegree(1,:,2),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    title('Average degree','Interpreter','latex')
    %if strcmp(simtype,'baseline')
       % legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %elseif strcmp(simtype,'crisis')
        %legend([h1 h2],{'FRFA OFF','FRFA ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %end
    xlabel('Iteration step','Interpreter','latex')
%----------------------------------------------------------------------------------------------------------    
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_F,'IBnetwork','.pdf'));

%-----------------------------------------------------
% Overlapping portfolio network
%-----------------------------------------------------

figure
subplot(1,3,1)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),10],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(mu_A(1,:,1),'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(mu_A(1,:,2),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    title('$\mu_{A}$ (Average diversification)','Interpreter','latex')
    if strcmp(simtype,'baseline')
        legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},'Location','best','FontSize',8,'Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        legend([h1 h2],{'FRFA OFF','FRFA ON'},'Location','best','FontSize',8,'Interpreter','latex')
    end
    xlabel('Iteration step','Interpreter','latex')
%----------------------------------------------------------------------------------------------------------    
subplot(1,3,2)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),10],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on;
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(mu_B(1,:,1),'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(mu_B(1,:,2),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    title('$\mu_{B}$ (Average diversification)','Interpreter','latex')
    %if strcmp(simtype,'baseline')
        %legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %elseif strcmp(simtype,'crisis')
        %legend([h1 h2],{'FRFA OFF','FRFA ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %end
    xlabel('Iteration step','Interpreter','latex')
%----------------------------------------------------------------------------------------------------------    
subplot(1,3,3)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),1],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(N_vec(1,:,1)./M_vec(1,:,1),'Color','r','LineWidth',1.1);
    hold on
    h2 = plot(N_vec(1,:,2)./M_vec(1,:,2),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    title('$N/M$ (Crowding)','Interpreter','latex')
    %if strcmp(simtype,'baseline')
        %legend([h1 h2],{'CB refinancing OFF','CB refinancing ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %elseif strcmp(simtype,'crisis')
        %legend([h1 h2],{'FRFA OFF','FRFA ON'},'Location','best','FontSize',8,'Interpreter','latex')
    %end
    xlabel('Iteration step','Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_F,'OPnetwork','.pdf'));

%-----------------------------------------------------
% Active banks/assets over simulation
%-----------------------------------------------------

figure
subplot(1,2,1)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),n_banks],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        %plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        %hold on;
    end
    h1 = plot(N_vec(1,:,1),'Color','r','LineWidth',1.1);
    hold on;
    h2 = plot(M_vec(1,:,1),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    if strcmp(simtype,'baseline')
        title('(a) CB refinancing OFF','Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title('(a) FRFA OFF','Interpreter','latex')
    end
    legend([h1 h2],{'\# of active banks','\# of active assets'},'Location','best','FontSize',8,'Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
subplot(1,2,2)
    if strcmp(simtype,'crisis')
        rectangle('Position',[shocktime(1),0,shocktime(end)-shocktime(1),n_banks],...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0.8 0.8 0.8]);
        hold on
        plot([FRFAtime(1), FRFAtime(1)],get(gca,'ylim'),'--k','LineWidth',2);
        hold on;
    end
    h1 = plot(N_vec(1,:,2),'Color','r','LineWidth',1.1);
    hold on;
    h2 = plot(M_vec(1,:,2),'Color','b','LineWidth',1.1);
    xlim([0 T_sim]);
    grid on;
    if strcmp(simtype,'baseline')
        title('(a) CB refinancing ON','Interpreter','latex')
    elseif strcmp(simtype,'crisis')
        title('(a) FRFA ON','Interpreter','latex')
    end
    %legend([h1 h2],{'\# of active banks','\# of active assets'},'Location','best','FontSize',8,'Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
%----------------------------------------------------------------------------------------------------------
set(gcf,'renderer','painters');
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3),pos(4)])
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 16, 8],'PaperUnits','Centimeters','PaperSize',[16,8])
print(gcf,'-dpdf',strcat(fig_output_F,'ActiveBanks+Assets','.pdf'));

end