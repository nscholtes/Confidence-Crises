function Results_balancesheet(fig_output,TOT_T_matrices,TOT_NT_matrices,normalizegraph)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Balance sheet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_output_BS = strcat(fig_output,'Results/');
 % Select whether to normalize results by dividing by total assets

if strcmp(normalizegraph,'Y')
    normfactor = TOT_T_matrices(1,:);
    normlabel  = 'norm_';
elseif strcmp(normalizegraph,'N')
    normfactor = ones(1,size(TOT_T_matrices,2));
    normlabel  = 'raw_';
end

% Aggregate balance sheet size

figure
plot(TOT_T_matrices(1,:),'LineWidth',1.1)
grid on;
xlabel('Iteration step','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_BS,'totalassets.pdf'));

% Evolution of balance sheet componenents
figure
subplot(2,2,1)
    plot(TOT_T_matrices(2,:)./normfactor,'LineWidth',1.1)
    grid on;
    title('Cash','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
subplot(2,2,2)
    plot(TOT_T_matrices(3,:)./normfactor,'LineWidth',1.1)
    grid on;
    title('External assets','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
subplot(2,2,3)
    plot(TOT_T_matrices(6,:)./normfactor,'LineWidth',1.1)
    grid on;
    title('Deposits','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')
subplot(2,2,4)
    plot(TOT_T_matrices(7,:)./normfactor,'LineWidth',1.1)
    grid on;
    title('Capital','Interpreter','latex')
    xlabel('Iteration step','Interpreter','latex')  
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_BS,normlabel,'balancesheet.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Investment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(TOT_T_matrices(4,:)./normfactor,'Color','b','LineWidth',1.1)
hold on
plot(TOT_T_matrices(5,:)./normfactor,'Color','r','LineWidth',1.1)
grid on;
%grid minor;
legend({'Total planned investment','Final investment'},'Location','best','FontSize',8,'Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_BS,normlabel,'investment.pdf'));


end
