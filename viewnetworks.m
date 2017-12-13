function[Networkdata,Graphdata] = viewnetworks(fig_output,networkparameters,a,timeslices_IBN,timeslices_OPN,UU_adjmats,DW_adjmats,BP_mats,tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_output_VN = strcat(fig_output,'Network/');

% Presetting loop parameters for subplots associated to selected timeslices
for i = 1:numel(timeslices_IBN)
    Slice_IBN{i} = strcat('Slice',num2str(timeslices_IBN(i)));
    graphtitle_IBN{i} = strcat('t = ',num2str(timeslices_IBN(i)));
end

for i = 1:numel(timeslices_OPN)
    Slice_OPN{i} = strcat('Slice',num2str(timeslices_OPN(i)));
    graphtitle_OPN{i} = strcat('t = ',num2str(timeslices_OPN(i)));
end

adj_mat   = UU_adjmats;

diadj_mat = DW_adjmats(:,:,:,1);
loan_mat  = DW_adjmats(:,:,:,2);

opnet_mat = BP_mats;

markernorm  = networkparameters(1);
edgewnorm   = networkparameters(2);  
graphlayout = networkparameters(3); 

for i = 1:numel(Slice_IBN)
    % Bank nodes
    Networkdata_nodes.banks.(Slice_IBN{i}).sizes_zero   = (a(:,i)./(max(a(:,1))/markernorm)); % Use largest bank at start of simulation to size nodes
    Networkdata_nodes.banks.(Slice_IBN{i}).nonzero_ids  = find(Networkdata_nodes.banks.(Slice_IBN{i}).sizes_zero~=0);
    Networkdata_nodes.banks.(Slice_IBN{i}).sizes_nozero = Networkdata_nodes.banks.(Slice_IBN{i}).sizes_zero;
    Networkdata_nodes.banks.(Slice_IBN{i}).sizes_nozero(Networkdata_nodes.banks.(Slice_IBN{i}).sizes_nozero==0) = [];
    Networkdata_nodes.banks.(Slice_IBN{i}).numnodes     = numel(Networkdata_nodes.banks.(Slice_IBN{i}).sizes_nozero);
    Networkdata_nodes.banks.(Slice_IBN{i}).nodearray    = 1:Networkdata_nodes.banks.(Slice_IBN{i}).numnodes;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Undirected, unweighted networks (UUN): Existing relationships between banks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numel(Slice_IBN)
    Networkdata_UU.(Slice_IBN{i}).adj_mat = adj_mat(find(Networkdata_nodes.banks.(Slice_IBN{i}).sizes_zero),find(Networkdata_nodes.banks.(Slice_IBN{i}).sizes_zero),i);
end

% Initial simulated network
G_UU.(Slice_IBN{1}) = graph(Networkdata_UU.(Slice_IBN{1}).adj_mat);
G_UU.(Slice_IBN{1}).Nodes.Size = Networkdata_nodes.banks.(Slice_IBN{1}).sizes_nozero;


figure
subplot(1,2,1)
    h = plot(G_UU.(Slice_IBN{1}),'Layout','force','MarkerSize',Networkdata_nodes.banks.(Slice_IBN{1}).sizes_nozero);
    labelnode(h,Networkdata_nodes.banks.(Slice_IBN{1}).nodearray,Networkdata_nodes.banks.(Slice_IBN{1}).nonzero_ids)
    set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
    title('(a) Force layout','Interpreter','latex')
subplot(1,2,2)
    h = plot(G_UU.(Slice_IBN{1}),'Layout',graphlayout,'MarkerSize',Networkdata_nodes.banks.(Slice_IBN{1}).sizes_nozero);
    labelnode(h,Networkdata_nodes.banks.(Slice_IBN{1}).nodearray,Networkdata_nodes.banks.(Slice_IBN{1}).nonzero_ids)
    set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
    title('(b) Circular layout','Interpreter','latex')
set(gcf,'renderer','painters');
hfig = tightfig;
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig,'-dpdf',strcat(fig_output_VN,'interbanknetwork_initial.pdf'));

% Evolving network structure at four time points over the simulation  figure
figure
for i = 2:numel(Slice_IBN)
    subplot(2,2,i-1)
        G_UU.(Slice_IBN{i}) = graph(Networkdata_UU.(Slice_IBN{i}).adj_mat);
        G_UU.(Slice_IBN{i}).Nodes.Size = Networkdata_nodes.banks.(Slice_IBN{i}).sizes_nozero;
        h = plot(G_UU.(Slice_IBN{i}),'Layout',graphlayout,'MarkerSize',Networkdata_nodes.banks.(Slice_IBN{i}).sizes_nozero);
        labelnode(h,Networkdata_nodes.banks.(Slice_IBN{i}).nodearray,Networkdata_nodes.banks.(Slice_IBN{i}).nonzero_ids)
        set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
        title(graphtitle_IBN{i},'Interpreter','latex')
end
set(gcf,'renderer','painters');
hfig = tightfig;
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig,'-dpdf',strcat(fig_output_VN,'interbanknetwork_evolution.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Directed, Weighted networks: Interbank market dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

for i = 1:numel(Slice_IBN)
    Networkdata_DW.(Slice_IBN{i}).diadj_mat = diadj_mat(:,:,i);  
    Networkdata_DW.(Slice_IBN{i}).loan_mat  = loan_mat(:,:,i);
    Networkdata_DW.(Slice_IBN{i}).edgelist = adj2edgeL(diadj_mat(:,:,i));
    if ~isempty(Networkdata_DW.(Slice_IBN{i}).edgelist)
        Networkdata_DW.(Slice_IBN{i}).srcnodes = Networkdata_DW.(Slice_IBN{i}).edgelist(:,1);
        Networkdata_DW.(Slice_IBN{i}).ternodes = Networkdata_DW.(Slice_IBN{i}).edgelist(:,2);
        Networkdata_DW.(Slice_IBN{i}).activenodes = unique([Networkdata_DW.(Slice_IBN{i}).srcnodes; Networkdata_DW.(Slice_IBN{i}).ternodes]);
        Networkdata_DW.(Slice_IBN{i}).numactivenodes = numel(Networkdata_DW.(Slice_IBN{i}).activenodes);
        Networkdata_DW.(Slice_IBN{i}).activenodearray = 1:numel(Networkdata_DW.(Slice_IBN{i}).activenodes);
        Networkdata_DW.(Slice_IBN{i}).idx         = sub2ind(size(Networkdata_DW.(Slice_IBN{i}).diadj_mat),...
        Networkdata_DW.(Slice_IBN{i}).srcnodes,Networkdata_DW.(Slice_IBN{i}).ternodes);
        Networkdata_DW.(Slice_IBN{i}).weights     = Networkdata_DW.(Slice_IBN{i}).loan_mat(Networkdata_DW.(Slice_IBN{i}).idx);
        Networkdata_DW.(Slice_IBN{i}).edgewidths = edgewnorm*Networkdata_DW.(Slice_IBN{i}).weights/max(Networkdata_DW.(Slice_IBN{i}).weights);
    end
end

Networkdata.nodes = Networkdata_nodes.banks;
Networkdata.UU    = Networkdata_UU;
Networkdata.DW    = Networkdata_DW;

% 1st period interbank dynamics
G_DW.(Slice_IBN{1}) = digraph(Networkdata_DW.(Slice_IBN{1}).srcnodes,Networkdata_DW.(Slice_IBN{1}).ternodes,Networkdata_DW.(Slice_IBN{1}).weights);
G_DW.(Slice_IBN{1}) = rmnode(G_DW.(Slice_IBN{1}),setdiff(Networkdata_nodes.banks.(Slice_IBN{1}).nodearray,Networkdata_DW.(Slice_IBN{1}).activenodes));
G_DW.(Slice_IBN{1}).Nodes.Size = Networkdata_nodes.banks.(Slice_IBN{1}).sizes_zero(Networkdata_DW.(Slice_IBN{1}).activenodes);

figure
subplot(1,2,1)
    h = plot(G_DW.(Slice_IBN{1}),'Layout','force','LineWidth',Networkdata_DW.(Slice_IBN{1}).edgewidths,'MarkerSize',G_DW.(Slice_IBN{1}).Nodes.Size);
    labelnode(h,Networkdata_DW.(Slice_IBN{1}).activenodearray,Networkdata_DW.(Slice_IBN{1}).activenodes)
    set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
    title('(a) Force layout','Interpreter','latex')
subplot(1,2,2)
    h = plot(G_DW.(Slice_IBN{1}),'Layout',graphlayout,'LineWidth',Networkdata_DW.(Slice_IBN{1}).edgewidths,'MarkerSize',G_DW.(Slice_IBN{1}).Nodes.Size);
    labelnode(h,Networkdata_DW.(Slice_IBN{1}).activenodearray,Networkdata_DW.(Slice_IBN{1}).activenodes)
    set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
    title('(b) Circular layout','Interpreter','latex')
set(gcf,'renderer','painters');
hfig = tightfig;
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig,'-dpdf',strcat(fig_output_VN,'interbankloans_initial.pdf'));

% Interbank dynamics for selected $timeslices$
figure
for i = 2:numel(Slice_IBN)
    subplot(2,numel(Slice_IBN)-3,i-1)
    if ~isempty(Networkdata_DW.(Slice_IBN{i}).edgelist)
        G_DW.(Slice_IBN{i}) = digraph(Networkdata_DW.(Slice_IBN{i}).srcnodes,Networkdata_DW.(Slice_IBN{i}).ternodes,Networkdata_DW.(Slice_IBN{i}).weights);
        G_DW.(Slice_IBN{i}) = rmnode(G_DW.(Slice_IBN{i}),setdiff(Networkdata_nodes.banks.(Slice_IBN{1}).nodearray,Networkdata_DW.(Slice_IBN{i}).activenodes));
        G_DW.(Slice_IBN{i}).Nodes.Size = Networkdata_nodes.banks.(Slice_IBN{i}).sizes_zero(Networkdata_DW.(Slice_IBN{i}).activenodes);
       % Networkdata_DW.(Slice{i}).activenodes
        %G_DW.(Slice{i}).Nodes.Size
        h = plot(G_DW.(Slice_IBN{i}),'Layout',graphlayout,'LineWidth',Networkdata_DW.(Slice_IBN{i}).edgewidths,'MarkerSize',G_DW.(Slice_IBN{i}).Nodes.Size+tol);
        labelnode(h,Networkdata_DW.(Slice_IBN{i}).activenodearray,Networkdata_DW.(Slice_IBN{i}).activenodes)
    else
        h = plot(G_UU.(Slice_IBN{i}),'Layout',graphlayout,'MarkerSize',Networkdata_nodes.banks.(Slice_IBN{i}).sizes_nozero);
        labelnode(h,Networkdata_nodes.banks.(Slice_IBN{i}).nodearray,Networkdata_nodes.banks.(Slice_IBN{i}).nonzero_ids)
    end
        set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
        title(graphtitle_IBN{i},'Interpreter','latex')
end 
set(gcf,'renderer','painters');
hfig = tightfig;
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig,'-dpdf',strcat(fig_output_VN,'interbankloans_evolution.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Weighted bipartite graph: Overlapping portfolio dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for i = 1:numel(Slice_OPN) 
   
    Networkdata_BP.(Slice_OPN{i}).opnet_mat_zero   = opnet_mat(:,:,i);
    Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero = Networkdata_BP.(Slice_OPN{i}).opnet_mat_zero;
    
    Networkdata_BP.(Slice_OPN{i}).bankids  = find(sum(Networkdata_BP.(Slice_OPN{i}).opnet_mat_zero,2)~=0);
    Networkdata_BP.(Slice_OPN{i}).assetids = find(sum(Networkdata_BP.(Slice_OPN{i}).opnet_mat_zero,1)~=0);
     
    Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero(~any(Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero,2),:) = [];
    Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero(:,~any(Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero,1)) = [];
    
    Networkdata_BP.(Slice_OPN{i}).numbanks   = size(Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero,1);
    Networkdata_BP.(Slice_OPN{i}).bankarray  = 1:size(Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero,1);
    Networkdata_BP.(Slice_OPN{i}).numassets  = size(Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero,2);
    Networkdata_BP.(Slice_OPN{i}).assetarray = 1:size(Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero,2);

% Represent bipartite graph as an [N+M]x[N+M] adjacency matrix
    % Populate TL quadrant: ZERO MATRIX - Banks in bipartite graph are not connected - NxN
    Networkdata_BP.(Slice_OPN{i}).opvizmat(Networkdata_BP.(Slice_OPN{i}).bankarray,...
        Networkdata_BP.(Slice_OPN{i}).bankarray) = zeros(Networkdata_BP.(Slice_OPN{i}).numbanks);
    % Populate BR quadrant: ZERO MATRIX - Assets in bipartite graph are not connected - MxM
    Networkdata_BP.(Slice_OPN{i}).opvizmat(Networkdata_BP.(Slice_OPN{i}).numbanks+1:(Networkdata_BP.(Slice_OPN{i}).numbanks + Networkdata_BP.(Slice_OPN{i}).numassets),...
        Networkdata_BP.(Slice_OPN{i}).numbanks+1:(Networkdata_BP.(Slice_OPN{i}).numbanks+Networkdata_BP.(Slice_OPN{i}).numassets)) = ...
         zeros(Networkdata_BP.(Slice_OPN{i}).numassets);
    % Populate BL quadrant: Asset-Bank links - MxN
    Networkdata_BP.(Slice_OPN{i}).opvizmat(Networkdata_BP.(Slice_OPN{i}).numbanks+1:(Networkdata_BP.(Slice_OPN{i}).numbanks+Networkdata_BP.(Slice_OPN{i}).numassets),...
        Networkdata_BP.(Slice_OPN{i}).bankarray) = (Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero)';
    % Populate TR quadrant:  Bank-Asset links - NxM
    Networkdata_BP.(Slice_OPN{i}).opvizmat(Networkdata_BP.(Slice_OPN{i}).bankarray,...
        Networkdata_BP.(Slice_OPN{i}).numbanks+1:(Networkdata_BP.(Slice_OPN{i}).numbanks+Networkdata_BP.(Slice_OPN{i}).numassets)) =...
        Networkdata_BP.(Slice_OPN{i}).opnet_mat_nozero;
end

% Initial overlapping portfolio network

G_BP.(Slice_OPN{1}) = graph(Networkdata_BP.(Slice_OPN{1}).opvizmat);

figure
h=plot(G_BP.(Slice_OPN{1}));
% Assets represented in red
highlight(h,Networkdata_BP.(Slice_OPN{1}).numbanks+1:(Networkdata_BP.(Slice_OPN{1}).numbanks+Networkdata_BP.(Slice_OPN{1}).numassets),'NodeColor','r');
% Labelling active BANKS
labelnode(h,Networkdata_BP.(Slice_OPN{1}).bankarray,Networkdata_BP.(Slice_OPN{1}).bankids)
% Labelling active ASSETS
labelnode(h,Networkdata_BP.(Slice_OPN{1}).numbanks+1:Networkdata_BP.(Slice_OPN{1}).numbanks+Networkdata_BP.(Slice_OPN{1}).numassets,...
    Networkdata_BP.(Slice_OPN{1}).assetids)
% Reorganising to visualise as bipartite graph
%%% Banks on LHS
h.XData(1:Networkdata_BP.(Slice_OPN{1}).numbanks) = 1;
h.YData(1:Networkdata_BP.(Slice_OPN{1}).numbanks) = linspace(0,1,Networkdata_BP.(Slice_OPN{1}).numbanks);
%h.YData(1:Networkdata_BP.(Slice{1}).numbanks) = linspace(0.5-(Networkdata_BP.(Slice{1}).numbanks/(2*Networkdata_BP.(Slice{1}).numassets)),...
    %0.5+(Networkdata_BP.(Slice{1}).numbanks/(2*Networkdata_BP.(Slice{1}).numassets)),Networkdata_BP.(Slice{1}).numbanks);
%%% Assets on RHS
h.XData(Networkdata_BP.(Slice_OPN{1}).numbanks+1:(Networkdata_BP.(Slice_OPN{1}).numbanks+Networkdata_BP.(Slice_OPN{1}).numassets)) = 2;
h.YData(Networkdata_BP.(Slice_OPN{1}).numbanks+1:(Networkdata_BP.(Slice_OPN{1}).numbanks+Networkdata_BP.(Slice_OPN{1}).numassets))...
    = linspace(0,1,Networkdata_BP.(Slice_OPN{1}).numassets);
set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
set(gcf,'renderer','painters');
hfig = tightfig;
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig,'-dpdf',strcat(fig_output_VN,'overlappingportfolio_initial.pdf'));

% Overlapping portfolio network evolution
figure
for i = 2:numel(Slice_OPN)
    subplot(3,numel(Slice_OPN)-5,i-1)
        G_BP.(Slice_OPN{i}) = graph(Networkdata_BP.(Slice_OPN{i}).opvizmat);
        h=plot(G_BP.(Slice_OPN{i}),'NodeLabel',[]);
        highlight(h,Networkdata_BP.(Slice_OPN{i}).numbanks+1:(Networkdata_BP.(Slice_OPN{i}).numbanks+Networkdata_BP.(Slice_OPN{i}).numassets),'NodeColor','r');
        
        %labelnode(h,Networkdata_BP.(Slice{i}).bankarray,Networkdata_BP.(Slice{i}).bankids)
        %labelnode(h,Networkdata_BP.(Slice{i}).numbanks+1:Networkdata_BP.(Slice{i}).numbanks+Networkdata_BP.(Slice{i}).numassets,...
            %Networkdata_BP.(Slice{i}).assetids)
        
% Reorganising to visualise as bipartite graph
        h.XData(1:Networkdata_BP.(Slice_OPN{i}).numbanks) = 1;
        %h.YData(1:Networkdata_BP.(Slice{i}).numbanks) = linspace(0.5-(Networkdata_BP.(Slice{i}).numbanks/(2*Networkdata_BP.(Slice{i}).numassets)),...
            %0.5+(Networkdata_BP.(Slice{i}).numbanks/(2*Networkdata_BP.(Slice{i}).numassets)),Networkdata_BP.(Slice{i}).numbanks);
        h.YData(1:Networkdata_BP.(Slice_OPN{i}).numbanks) = linspace(0,1,Networkdata_BP.(Slice_OPN{i}).numbanks);
        h.XData(Networkdata_BP.(Slice_OPN{i}).numbanks+1:(Networkdata_BP.(Slice_OPN{i}).numbanks+Networkdata_BP.(Slice_OPN{i}).numassets)) = 2;
        h.YData(Networkdata_BP.(Slice_OPN{i}).numbanks+1:(Networkdata_BP.(Slice_OPN{i}).numbanks+Networkdata_BP.(Slice_OPN{i}).numassets))...
            = linspace(0,1,Networkdata_BP.(Slice_OPN{i}).numassets);
        title(graphtitle_OPN{i},'Interpreter','latex')
        set(gca,'XLim',[1 2],'YLim',[0 1])
        set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);    
end
set(gcf,'renderer','painters');
%hfig = tightfig;
%set(gcf,'Units','Inches');
%pos = get(hfig,'Position');
set(gcf,'Units', 'Centimeters', 'Position', [0, 0, 18, 23],'PaperUnits','Centimeters','PaperSize',[18,23])
%set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_VN,'overlappingportfolio_evolution.pdf'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Output collected network information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

Networkdata.nodes = Networkdata_nodes.banks;
Networkdata.UU    = Networkdata_UU;
Networkdata.DW    = Networkdata_DW;
Networkdata.BP    = Networkdata_BP;

Graphdata.UU = G_UU; 
Graphdata.DW = G_DW;
Graphdata.BP = G_BP;

end