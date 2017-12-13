function [EONIA_data,MRO_data,MFI_data] = dataoutput()

% Set directory for input data and figure output

data_dir = '/Users/nscholte/Desktop/Research/Ch.1 - Confidence crises/Data/';
fig_output_data = '/Users/nscholte/Desktop/Research/Ch.1 - Confidence crises/Drafts/Current version/Figures/Data/';

% Key dates

crisis_start_date = datenum('2008-09-15');
FRFA_start_date   = datenum('2008-10-15');

x_axis_datestring = {'2008-06-01','2008-07-01','2008-08-01','2008-09-01','2008-10-01','2008-11-01','2008-12-01','2009-01-01'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. EONIA volumes and rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Read in data
%--------------------------------------------------------------------------

% Volumes
fid = fopen(strcat(data_dir,'EONIA_volumes.csv'));
EONIA_volume_data = textscan(fid, '%s %f', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);
EONIA_volume_data{1}  = datenum(EONIA_volume_data{1});
fin_EONIA_volume_data = cell2mat(EONIA_volume_data);

% Rates
fid = fopen(strcat(data_dir,'EONIA_rates.csv'));
EONIA_rate_data = textscan(fid, '%s %f', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);
EONIA_rate_data{1}  = datenum(EONIA_rate_data{1});
fin_EONIA_rate_data = cell2mat(EONIA_rate_data);

% Store in same structure
EONIA_data         = struct;
EONIA_data.volumes = fin_EONIA_volume_data;
EONIA_data.rates   = fin_EONIA_rate_data;

%--------------------------------------------------------------------------
% Plot data
%--------------------------------------------------------------------------

figure
subplot(1,2,1)
    plot(fin_EONIA_volume_data(:,1),fin_EONIA_volume_data(:,2),'LineWidth',1.1)
    hold on
    plot([crisis_start_date,crisis_start_date],get(gca,'ylim'),'--k','LineWidth',2);
    hold on
    plot([FRFA_start_date,FRFA_start_date],get(gca,'ylim'),'--r','LineWidth',2);
    dateFormat = 19;
    ax = gca;
    ax.XTick = datenum(x_axis_datestring);
    xlim([fin_EONIA_volume_data(end,1),fin_EONIA_volume_data(1,1)])
    datetick('x',dateFormat,'keepticks')
    xtickangle(45)
    grid on
    ax.YAxis.Exponent = 0;
    ylabel('Volumes [EUR Million]','Interpreter','latex')
subplot(1,2,2)
    plot(fin_EONIA_rate_data(:,1),fin_EONIA_rate_data(:,2),'LineWidth',1.1)
    hold on
    plot([crisis_start_date,crisis_start_date],get(gca,'ylim'),'--k','LineWidth',2);
    hold on
    plot([FRFA_start_date,FRFA_start_date],get(gca,'ylim'),'--r','LineWidth',2);
    dateFormat = 19;
    ax = gca;
    ax.XTick = datenum(x_axis_datestring);
    xlim([fin_EONIA_rate_data(end,1),fin_EONIA_rate_data(1,1)])
    datetick('x',dateFormat,'keepticks')
    xtickangle(45)
    grid on
    ax.YAxis.Exponent = 0;
    ylabel('Rates [\%]','Interpreter','latex')    
set(gcf,'renderer','painters');
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 18, 8],'PaperUnits','Centimeters','PaperSize',[18,8])
print(gcf,'-dpdf',strcat(fig_output_data,'EONIA','.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. ECB allotment volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Read in data
%--------------------------------------------------------------------------

% MRO - VRT
fid = fopen(strcat(data_dir,'ECB_allotments_MRO_VRT.csv'));
ECB_MRO_VRT_data = textscan(fid, '%s %f %f %f %f %f %f %f', 'Delimiter', ',');
fclose(fid);
ECB_MRO_VRT_data{1}  = datenum(ECB_MRO_VRT_data{1});
fin_ECB_MRO_VRT_data = cell2mat(ECB_MRO_VRT_data);

% MRO - FRT
fid = fopen(strcat(data_dir,'ECB_allotments_MRO_FRT.csv'));
ECB_MRO_FRT_data = textscan(fid, '%s %f %f %f %f %f %f', 'Delimiter', ',');
fclose(fid);
ECB_MRO_FRT_data{1}  = datenum(ECB_MRO_FRT_data{1});
fin_ECB_MRO_FRT_data = cell2mat(ECB_MRO_FRT_data);

%LTRO
fid = fopen(strcat(data_dir,'ECB_allotments_LTRO.csv'));
ECB_LTRO_data = textscan(fid, '%s %f %f %f %f %f %f', 'Delimiter', ',');
fclose(fid);
ECB_LTRO_data{1}  = datenum(ECB_LTRO_data{1});
fin_ECB_LTRO_data = cell2mat(ECB_LTRO_data);


% Store in same structure
MRO_data      = struct;
MRO_data.VRT  = fin_ECB_MRO_VRT_data;
MRO_data.FRT  = fin_ECB_MRO_FRT_data;
MRO_data.LTRO = fin_ECB_LTRO_data;

%--------------------------------------------------------------------------
% Plot data
%--------------------------------------------------------------------------

% Main refinancing operations
figure
subplot(1,2,1)  
    % MRO - VRT: Bids
    h1 = plot(fin_ECB_MRO_VRT_data(:,1),fin_ECB_MRO_VRT_data(:,2),'-r','LineWidth',1.1);
    hold on;
    % MRO - VRT: Allotment
    h2 = plot(fin_ECB_MRO_VRT_data(:,1),fin_ECB_MRO_VRT_data(:,4),'-b','LineWidth',1.1);
    hold on;
    % MRO - FRT
    h3 = plot(fin_ECB_MRO_FRT_data(:,1),fin_ECB_MRO_FRT_data(:,2),'-g','LineWidth',1.1);
    hold on
    plot([crisis_start_date,crisis_start_date],get(gca,'ylim'),'--k','LineWidth',2);
    hold on
    plot([FRFA_start_date,FRFA_start_date],get(gca,'ylim'),'--r','LineWidth',2);
    dateFormat = 19;
    ax = gca;
    ax.XTick = datenum(x_axis_datestring);
    %xlim([fin_EONIA_volume_data(end,1),fin_EONIA_volume_data(1,1)])
    datetick('x',dateFormat,'keepticks')
    xtickangle(45)
    grid on
    ax.YAxis.Exponent = 0;
    legend([h1 h2 h3],{'VRT: Bid amount','VRT: Alloted amount','FRT: Alloted amount'},...
         'Location','best','FontSize',6,'Interpreter','latex')
    ylabel('MRO volumes [EUR million]','Interpreter','latex')
subplot(1,2,2)
    % MRO - VRT: Number of bids
    h1 = plot(fin_ECB_MRO_VRT_data(:,1),fin_ECB_MRO_VRT_data(:,3),'-b','LineWidth',1.1);
    hold on;
    % MRO - FRT: Number of bids
    h2 = plot(fin_ECB_MRO_FRT_data(:,1),fin_ECB_MRO_FRT_data(:,3),'-r','LineWidth',1.1);
    hold on;
    plot([crisis_start_date,crisis_start_date],get(gca,'ylim'),'--k','LineWidth',2);
    hold on
    plot([FRFA_start_date,FRFA_start_date],get(gca,'ylim'),'--r','LineWidth',2);
    dateFormat = 19;
    ax = gca;
    ax.XTick = datenum(x_axis_datestring);
    %xlim([fin_EONIA_volume_data(end,1),fin_EONIA_volume_data(1,1)])
    datetick('x',dateFormat,'keepticks')
    xtickangle(45)
    grid on
    ax.YAxis.Exponent = 0;
    legend([h1 h2],{'VRT','FRT'},...
                'Location','best','FontSize',6,'Interpreter','latex')
    ylabel('Number of bids','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 18, 8],'PaperUnits','Centimeters','PaperSize',[18,8])
print(gcf,'-dpdf',strcat(fig_output_data,'ECB_MRO','.pdf'));

% Long-Term Refinancing Operations
figure
subplot(1,2,1)
    % LTRO: Bids
    h1 = plot(fin_ECB_LTRO_data(:,1),fin_ECB_LTRO_data(:,2),'-or','LineWidth',1.1);
    hold on
    % LTRO: Allotment
    h2 = plot(fin_ECB_LTRO_data(:,1),fin_ECB_LTRO_data(:,4),'-b','LineWidth',1.1);
    hold on
    plot([crisis_start_date,crisis_start_date],get(gca,'ylim'),'--k','LineWidth',2);
    hold on
    plot([FRFA_start_date,FRFA_start_date],get(gca,'ylim'),'--r','LineWidth',2);
    dateFormat = 19;
    ax = gca;
    ax.XTick = datenum(x_axis_datestring);
    %xlim([fin_EONIA_volume_data(end,1),fin_EONIA_volume_data(1,1)])
    datetick('x',dateFormat,'keepticks')
    xtickangle(45)
    grid on
    ax.YAxis.Exponent = 0;
    legend([h1 h2],{'LTRO: Bids','LTRO: Allotment'},...
                'Location','best','FontSize',6,'Interpreter','latex')
    ylabel('LTRO volumes [EUR Million]','Interpreter','latex')
subplot(1,2,2)
    plot(fin_ECB_LTRO_data(:,1),fin_ECB_LTRO_data(:,3),'-r','LineWidth',1.1);
    hold on
    plot([crisis_start_date,crisis_start_date],get(gca,'ylim'),'--k','LineWidth',2);
    hold on
    plot([FRFA_start_date,FRFA_start_date],get(gca,'ylim'),'--r','LineWidth',2);
    dateFormat = 19;
    ax = gca;
    ax.XTick = datenum(x_axis_datestring);
    %xlim([fin_EONIA_volume_data(end,1),fin_EONIA_volume_data(1,1)])
    datetick('x',dateFormat,'keepticks')
    xtickangle(45)
    grid on
    ylabel('Number of bids','Interpreter','latex')
set(gcf,'renderer','painters');
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 18, 8],'PaperUnits','Centimeters','PaperSize',[18,8])
print(gcf,'-dpdf',strcat(fig_output_data,'ECB_LTRO','.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Balance sheet dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Read in data
%--------------------------------------------------------------------------

% MFI loans
fid = fopen(strcat(data_dir,'MFI_loans.csv'));
MFI_loans_data = textscan(fid, '%s %f', 'Delimiter', ',','HeaderLines', 1);
fclose(fid);
MFI_loans_data{1}  = datenum(MFI_loans_data{1});
fin_MFI_loans_data = cell2mat(MFI_loans_data);

% External assets
fid = fopen(strcat(data_dir,'MFI_Externalassets.csv'));
MFI_EA_data = textscan(fid, '%s %f', 'Delimiter', ',','HeaderLines', 1);
fclose(fid);
MFI_EA_data{1}     = datenum(MFI_EA_data{1});
fin_MFI_EA_data = cell2mat(MFI_EA_data);

MFI_data = struct;
MFI_data.loans = fin_MFI_loans_data;
MFI_data.EA = fin_MFI_EA_data;

%--------------------------------------------------------------------------
% Plot data
%--------------------------------------------------------------------------

x_axis_datestring = {'Q1-2008','Q2-2008','Q3-2008','Q4-2008','Q1-2009'};

figure
subplot(1,2,1)
    plot(fin_MFI_loans_data(:,1),fin_MFI_loans_data(:,2),'LineWidth',1.1)
    hold on
    plot([crisis_start_date,crisis_start_date],get(gca,'ylim'),'--k','LineWidth',2);
    hold on
    plot([FRFA_start_date,FRFA_start_date],get(gca,'ylim'),'--r','LineWidth',2);
    dateFormat = 19;
    ax = gca;
    ax.XTick = datenum(x_axis_datestring,'QQ-yyyy');
    %xlim([fin_EONIA_volume_data(end,1),fin_EONIA_volume_data(1,1)])
    datetick('x',dateFormat,'keepticks')
    xtickangle(45)
    grid on
subplot(1,2,2)
    plot(fin_MFI_EA_data(:,1),fin_MFI_EA_data(:,2),'LineWidth',1.1)
    hold on
    plot([crisis_start_date,crisis_start_date],get(gca,'ylim'),'--k','LineWidth',2);
    hold on
    plot([FRFA_start_date,FRFA_start_date],get(gca,'ylim'),'--r','LineWidth',2);
    dateFormat = 19;
    ax = gca;
    ax.XTick = datenum(x_axis_datestring,'QQ-yyyy');
    %xlim([fin_EONIA_volume_data(end,1),fin_EONIA_volume_data(1,1)])
    datetick('x',dateFormat,'keepticks')
    xtickangle(45)
    grid on
set(gcf,'renderer','painters');
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 18, 8],'PaperUnits','Centimeters','PaperSize',[18,8])
print(gcf,'-dpdf',strcat(fig_output_data,'MFI_BS','.pdf'));



end