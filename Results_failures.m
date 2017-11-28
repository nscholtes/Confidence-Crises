function [banksfail] = Results_failures(fig_output,banksfail,...
    FailCount,cum_Fails,capitalshortfall,N_vec,M_vec,numedges,density,avdegree,mu_A,mu_B,allsimlabels)

fig_output_F = strcat(fig_output,'Results/failures/');

FailCount_av  = FailCount(1,:); FailCount_min = FailCount(2,:); FailCount_max = FailCount(3,:);
cum_Fails_av  = cum_Fails(1,:); cum_Fails_min = cum_Fails(2,:); cum_Fails_max = cum_Fails(3,:);

capitalshortfall_av  = capitalshortfall(1,:); capitalshortfal_min = capitalshortfall(2,:); capitalshortfal_max = capitalshortfall(3,:);

N_av  = N_vec(1,:); N_min = N_vec(2,:); N_max = N_vec(3,:);
M_av  = M_vec(1,:); M_min = M_vec(2,:); M_max = M_vec(3,:);

numedges_av = numedges(1,:); numedges_min = numedges(2,:); numedges_max = numedges(3,:);
density_av  = density(1,:);  density_min  = density(2,:);  density_max  = density(3,:);
avdegree_av = avdegree(1,:); avdegree_min = avdegree(2,:); avdegree_max = avdegree(3,:);

mu_A_av = mu_A(1,:); mu_A_min = mu_A(2,:); mu_A_max = mu_A(3,:);
mu_B_av = mu_B(1,:); mu_B_min = mu_B(2,:); mu_B_max = mu_B(3,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Number of failures and capital shortfall over simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

figure
subplot(1,2,1)
    plot(FailCount_av)
    hold on
    plot(cum_Fails_av)
    xlabel('Iteration step','Interpreter','latex')
    grid on;
    legend({'\# of Failures','\# of Failures (cumulative)'},'Location','best','FontSize',8,'Interpreter','latex')
subplot(1,2,2)
    plot(capitalshortfall_av)
    hold on
    plot(cumsum(capitalshortfall_av))
    xlabel('Iteration step','Interpreter','latex')
    grid on;
    legend({'Capital shortfall','Capital shortfall (cumulative)'},'Location','best','FontSize',8,'Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_F,'failures','_',allsimlabels,'.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Graphing changing network structure due to ABM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

figure
subplot(1,3,1)
    plot(numedges_av)
    grid on;
    title('Number of edges','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
subplot(1,3,2)
    plot(density_av)
    grid on;
    title('Network density','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
subplot(1,3,3)
    plot(avdegree_av)
    grid on;
    title('Average degree','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_F,'IBnetwork','_',allsimlabels,'.pdf'));

figure
subplot(2,2,1)
    plot(mu_A_av)
    grid on;
    title('$\mu_{A}$ (Average diversification)','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
subplot(2,2,2)
    plot(mu_B_av)
    grid on;
    title('$\mu_{B}$','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
subplot(2,2,3)
    plot(N_av)
    hold on;
    plot(M_av)
    grid on;
    legend({'\# of active banks','\# of active assets'},'Location','best','FontSize',8,'Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
subplot(2,2,4)
    plot(N_av./M_av)
    grid on;
    title('$N/M$ (Crowding parameter)','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_F,'OPnetwork','_',allsimlabels,'.pdf'));

end