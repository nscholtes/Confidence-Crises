function Results_IBM(fig_output,TOT_T_matrices,TOT_NT_matrices,IBratemat,T)

fig_output_IBM = strcat(fig_output,'Results/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interbank market %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Volumes

figure
plot(TOT_T_matrices(8,:),'Color','b','LineWidth',1.1)
hold on
plot(TOT_T_matrices(9,:),'Color','r','LineWidth',1.1)
hold on
plot(TOT_T_matrices(10,:),'Color','k','LineWidth',1.1)
grid on;
xlabel('Iteration step','Interpreter','latex')
legend({'Requests','Hoarding','Loans'},'Location','best','FontSize',8)
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBM_phase1.pdf'));

figure
subplot(1,2,1)
    plot(TOT_T_matrices(11,:),'Color','b','LineWidth',1.1)
    hold on
    plot(TOT_T_matrices(12,:),'Color','r','LineWidth',1.1)
    grid on;
    legend({'Principal + interest','Repayment'},'Location','best','FontSize',8)
    xlabel('Iteration step','Interpreter','latex')
    title('(a)','Interpreter','latex')
subplot(1,2,2)
    plot(TOT_T_matrices(11,:)-TOT_T_matrices(12,:),'Color','k','LineWidth',1.1)
    xlabel('Iteration step','Interpreter','latex')
    title('(b)','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBM_phase2.pdf'));

% Rates

for t =1:T
    meanIBrate(t) = mean(nonzeros(IBratemat(:,:,t)));
    %minIBrate(t) =  min(nonzeros(IBratemat(:,:,t)));
    %maxIBrate(t) =  max(nonzeros(IBratemat(:,:,t)));    
end

figure
plot(meanIBrate,'Color','k','LineWidth',1.1)
%hold on
%plot(minIBrate,'Color','r','LineStyle','--','LineWidth',1.1);
%hold on
%plot(maxIBrate,'Color','b','LineStyle','-.','LineWidth',1.1);
xlabel('Iteration step','Interpreter','latex')
grid on;
%legend({'$\bar{r}^{b,t}$','$min\{r^{b,t}\}$','$max\{r^{b,t}\}$'},'Location','best','FontSize',8,'Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBM,'IBrates.pdf'));

end