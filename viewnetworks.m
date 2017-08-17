function[Networkdata,Graphdata] = viewnetworks(fig_output,networkparameters,a,timeslices,UU_adjmats,DW_adjmats,BP_mats)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_output_VN = strcat(fig_output,'Network/');

% Presetting loop parameters for subplots associated to selected timeslices
for i = 1:numel(timeslices)
    Slice{i} = strcat('Slice',num2str(timeslices(i)));
    graphtitle{i} = strcat('t = ',num2str(timeslices(i)));
end

adj_mat   = UU_adjmats;

diadj_mat = DW_adjmats(:,:,:,1);
loan_mat  = DW_adjmats(:,:,:,2);

opnet_mat = BP_mats;

markernorm  = networkparameters(1);
edgewnorm   = networkparameters(2);  
graphlayout = networkparameters(3); 

for i = 1:numel(Slice)
    % Bank nodes
    Networkdata_nodes.banks.(Slice{i}).sizes_zero   = (a(:,i)./(max(a(:,i))/markernorm)); % Use largest bank at start of simulation to size nodes
    Networkdata_nodes.banks.(Slice{i}).nonzero_ids  = find(Networkdata_nodes.banks.(Slice{i}).sizes_zero~=0);
    Networkdata_nodes.banks.(Slice{i}).sizes_nozero = Networkdata_nodes.banks.(Slice{i}).sizes_zero;
    Networkdata_nodes.banks.(Slice{i}).sizes_nozero(Networkdata_nodes.banks.(Slice{i}).sizes_nozero==0) = [];
    Networkdata_nodes.banks.(Slice{i}).numnodes     = numel(Networkdata_nodes.banks.(Slice{i}).sizes_nozero);
    Networkdata_nodes.banks.(Slice{i}).nodearray    = 1:Networkdata_nodes.banks.(Slice{i}).numnodes;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Undirected, unweighted networks (UUN): Existing relationships between banks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numel(Slice)
    Networkdata_UU.(Slice{i}).adj_mat = adj_mat(find(Networkdata_nodes.banks.(Slice{i}).sizes_zero),find(Networkdata_nodes.banks.(Slice{i}).sizes_zero),i);
end

% Initial simulated network
G_UU.(Slice{1}) = graph(Networkdata_UU.(Slice{1}).adj_mat);
G_UU.(Slice{1}).Nodes.Size = Networkdata_nodes.banks.(Slice{1}).sizes_nozero;


figure
subplot(1,2,1)
    h = plot(G_UU.(Slice{1}),'Layout','force','MarkerSize',Networkdata_nodes.banks.(Slice{1}).sizes_nozero);
    labelnode(h,Networkdata_nodes.banks.(Slice{1}).nodearray,Networkdata_nodes.banks.(Slice{1}).nonzero_ids)
    set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
    title('(a) Force layout','Interpreter','latex')
subplot(1,2,2)
    h = plot(G_UU.(Slice{1}),'Layout',graphlayout,'MarkerSize',Networkdata_nodes.banks.(Slice{1}).sizes_nozero);
    labelnode(h,Networkdata_nodes.banks.(Slice{1}).nodearray,Networkdata_nodes.banks.(Slice{1}).nonzero_ids)
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
for i = 2:numel(Slice)
    subplot(1,numel(Slice)-1,i-1)
        G_UU.(Slice{i}) = graph(Networkdata_UU.(Slice{i}).adj_mat);
        G_UU.(Slice{i}).Nodes.Size = Networkdata_nodes.banks.(Slice{i}).sizes_nozero;
        h = plot(G_UU.(Slice{i}),'Layout',graphlayout,'MarkerSize',Networkdata_nodes.banks.(Slice{i}).sizes_nozero);
        labelnode(h,Networkdata_nodes.banks.(Slice{i}).nodearray,Networkdata_nodes.banks.(Slice{i}).nonzero_ids)
        set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
        title(graphtitle{i},'Interpreter','latex')
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

for i = 1:numel(Slice)
    Networkdata_DW.(Slice{i}).diadj_mat = diadj_mat(:,:,i);  
    Networkdata_DW.(Slice{i}).loan_mat  = loan_mat(:,:,i);
    Networkdata_DW.(Slice{i}).edgelist = adj2edgeL(diadj_mat(:,:,i));
    if ~isempty(Networkdata_DW.(Slice{i}).edgelist)
        Networkdata_DW.(Slice{i}).srcnodes = Networkdata_DW.(Slice{i}).edgelist(:,1);
        Networkdata_DW.(Slice{i}).ternodes = Networkdata_DW.(Slice{i}).edgelist(:,2);
        Networkdata_DW.(Slice{i}).activenodes = unique([Networkdata_DW.(Slice{i}).srcnodes; Networkdata_DW.(Slice{i}).ternodes]);
        Networkdata_DW.(Slice{i}).numactivenodes = numel(Networkdata_DW.(Slice{i}).activenodes);
        Networkdata_DW.(Slice{i}).activenodearray = 1:numel(Networkdata_DW.(Slice{i}).activenodes);
        Networkdata_DW.(Slice{i}).idx         = sub2ind(size(Networkdata_DW.(Slice{i}).diadj_mat),...
        Networkdata_DW.(Slice{i}).srcnodes,Networkdata_DW.(Slice{i}).ternodes);
        Networkdata_DW.(Slice{i}).weights     = Networkdata_DW.(Slice{i}).loan_mat(Networkdata_DW.(Slice{i}).idx);
        Networkdata_DW.(Slice{i}).edgewidths = edgewnorm*Networkdata_DW.(Slice{i}).weights/max(Networkdata_DW.(Slice{i}).weights);
    end
end

Networkdata.nodes = Networkdata_nodes.banks;
Networkdata.UU    = Networkdata_UU;
Networkdata.DW    = Networkdata_DW;

% 1st period interbank dynamics
G_DW.(Slice{1}) = digraph(Networkdata_DW.(Slice{1}).srcnodes,Networkdata_DW.(Slice{1}).ternodes,Networkdata_DW.(Slice{1}).weights);
G_DW.(Slice{1}) = rmnode(G_DW.(Slice{1}),setdiff(Networkdata_nodes.banks.(Slice{1}).nodearray,Networkdata_DW.(Slice{1}).activenodes));
G_DW.(Slice{1}).Nodes.Size = Networkdata_nodes.banks.(Slice{1}).sizes_zero(Networkdata_DW.(Slice{1}).activenodes);

figure
subplot(1,2,1)
    h = plot(G_DW.(Slice{1}),'Layout','force','LineWidth',Networkdata_DW.(Slice{1}).edgewidths,'MarkerSize',G_DW.(Slice{1}).Nodes.Size);
    labelnode(h,Networkdata_DW.(Slice{1}).activenodearray,Networkdata_DW.(Slice{1}).activenodes)
    set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
    title('(a) Force layout','Interpreter','latex')
subplot(1,2,2)
    h = plot(G_DW.(Slice{1}),'Layout',graphlayout,'LineWidth',Networkdata_DW.(Slice{1}).edgewidths,'MarkerSize',G_DW.(Slice{1}).Nodes.Size);
    labelnode(h,Networkdata_DW.(Slice{1}).activenodearray,Networkdata_DW.(Slice{1}).activenodes)
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
for i = 2:numel(Slice)
    subplot(1,numel(Slice)-1,i-1)
    if ~isempty(Networkdata_DW.(Slice{i}).edgelist)
        G_DW.(Slice{i}) = digraph(Networkdata_DW.(Slice{i}).srcnodes,Networkdata_DW.(Slice{i}).ternodes,Networkdata_DW.(Slice{i}).weights);
        G_DW.(Slice{i}) = rmnode(G_DW.(Slice{i}),setdiff(Networkdata_nodes.banks.(Slice{1}).nodearray,Networkdata_DW.(Slice{i}).activenodes));
        G_DW.(Slice{i}).Nodes.Size = Networkdata_nodes.banks.(Slice{i}).sizes_zero(Networkdata_DW.(Slice{i}).activenodes);
        h = plot(G_DW.(Slice{i}),'Layout',graphlayout,'LineWidth',Networkdata_DW.(Slice{i}).edgewidths,'MarkerSize',G_DW.(Slice{i}).Nodes.Size);
        labelnode(h,Networkdata_DW.(Slice{i}).activenodearray,Networkdata_DW.(Slice{i}).activenodes)
    else
        h = plot(G_UU.(Slice{i}),'Layout',graphlayout,'MarkerSize',Networkdata_nodes.banks.(Slice{i}).sizes_nozero);
        labelnode(h,Networkdata_nodes.banks.(Slice{i}).nodearray,Networkdata_nodes.banks.(Slice{i}).nonzero_ids)
    end
        set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
        title(graphtitle{i},'Interpreter','latex')
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

for i = 1:numel(Slice)
    
    Networkdata_BP.(Slice{i}).opnet_mat_zero   = opnet_mat(:,:,i);
    Networkdata_BP.(Slice{i}).opnet_mat_nozero = Networkdata_BP.(Slice{i}).opnet_mat_zero;
    
    Networkdata_BP.(Slice{i}).bankids  = find(sum(Networkdata_BP.(Slice{i}).opnet_mat_zero,2)~=0);
    Networkdata_BP.(Slice{i}).assetids = find(sum(Networkdata_BP.(Slice{i}).opnet_mat_zero,1)~=0);
     
    Networkdata_BP.(Slice{i}).opnet_mat_nozero(~any(Networkdata_BP.(Slice{i}).opnet_mat_nozero,2),:) = [];
    Networkdata_BP.(Slice{i}).opnet_mat_nozero(:,~any(Networkdata_BP.(Slice{i}).opnet_mat_nozero,1)) = [];
    
    Networkdata_BP.(Slice{i}).numbanks   = size(Networkdata_BP.(Slice{i}).opnet_mat_nozero,1);
    Networkdata_BP.(Slice{i}).bankarray  = 1:size(Networkdata_BP.(Slice{i}).opnet_mat_nozero,1);
    Networkdata_BP.(Slice{i}).numassets  = size(Networkdata_BP.(Slice{i}).opnet_mat_nozero,2);
    Networkdata_BP.(Slice{i}).assetarray = 1:size(Networkdata_BP.(Slice{i}).opnet_mat_nozero,2);

% Represent bipartite graph as an [N+M]x[N+M] adjacency matrix
    % Populate TL quadrant: ZERO MATRIX - Banks in bipartite graph are not connected - NxN
    Networkdata_BP.(Slice{i}).opvizmat(Networkdata_BP.(Slice{i}).bankarray,...
        Networkdata_BP.(Slice{i}).bankarray) = zeros(Networkdata_BP.(Slice{i}).numbanks);
    % Populate BR quadrant: ZERO MATRIX - Assets in bipartite graph are not connected - MxM
    Networkdata_BP.(Slice{i}).opvizmat(Networkdata_BP.(Slice{i}).numbanks+1:(Networkdata_BP.(Slice{i}).numbanks + Networkdata_BP.(Slice{i}).numassets),...
        Networkdata_BP.(Slice{i}).numbanks+1:(Networkdata_BP.(Slice{i}).numbanks+Networkdata_BP.(Slice{i}).numassets)) = ...
         zeros(Networkdata_BP.(Slice{i}).numassets);
    % Populate BL quadrant: Asset-Bank links - MxN
    Networkdata_BP.(Slice{i}).opvizmat(Networkdata_BP.(Slice{i}).numbanks+1:(Networkdata_BP.(Slice{i}).numbanks+Networkdata_BP.(Slice{i}).numassets),...
        Networkdata_BP.(Slice{i}).bankarray) = (Networkdata_BP.(Slice{i}).opnet_mat_nozero)';
    % Populate TR quadrant:  Bank-Asset links - NxM
    Networkdata_BP.(Slice{i}).opvizmat(Networkdata_BP.(Slice{i}).bankarray,...
        Networkdata_BP.(Slice{i}).numbanks+1:(Networkdata_BP.(Slice{i}).numbanks+Networkdata_BP.(Slice{i}).numassets)) =...
        Networkdata_BP.(Slice{i}).opnet_mat_nozero;
end

% Initial overlapping portfolio network

G_BP.(Slice{1}) = graph(Networkdata_BP.(Slice{1}).opvizmat);

figure
h=plot(G_BP.(Slice{1}));
% Assets represented in red
highlight(h,Networkdata_BP.(Slice{1}).numbanks+1:(Networkdata_BP.(Slice{1}).numbanks+Networkdata_BP.(Slice{1}).numassets),'NodeColor','r');
% Labelling active BANKS
labelnode(h,Networkdata_BP.(Slice{1}).bankarray,Networkdata_BP.(Slice{1}).bankids)
% Labelling active ASSETS
labelnode(h,Networkdata_BP.(Slice{1}).numbanks+1:Networkdata_BP.(Slice{1}).numbanks+Networkdata_BP.(Slice{1}).numassets,...
    Networkdata_BP.(Slice{1}).assetids)
% Reorganising to visualise as bipartite graph
%%% Banks on LHS
h.XData(1:Networkdata_BP.(Slice{1}).numbanks) = 1;
h.YData(1:Networkdata_BP.(Slice{1}).numbanks) = linspace(0,1,Networkdata_BP.(Slice{1}).numbanks);
%h.YData(1:Networkdata_BP.(Slice{1}).numbanks) = linspace(0.5-(Networkdata_BP.(Slice{1}).numbanks/(2*Networkdata_BP.(Slice{1}).numassets)),...
    %0.5+(Networkdata_BP.(Slice{1}).numbanks/(2*Networkdata_BP.(Slice{1}).numassets)),Networkdata_BP.(Slice{1}).numbanks);
%%% Assets on RHS
h.XData(Networkdata_BP.(Slice{1}).numbanks+1:(Networkdata_BP.(Slice{1}).numbanks+Networkdata_BP.(Slice{1}).numassets)) = 2;
h.YData(Networkdata_BP.(Slice{1}).numbanks+1:(Networkdata_BP.(Slice{1}).numbanks+Networkdata_BP.(Slice{1}).numassets))...
    = linspace(0,1,Networkdata_BP.(Slice{1}).numassets);
set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);
set(gcf,'renderer','painters');
hfig = tightfig;
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig,'-dpdf',strcat(fig_output_VN,'overlappingportfolio_initial.pdf'));

% Overlapping portfolio network evolution
figure
for i = 2:numel(Slice)
    subplot(2,numel(Slice)-3,i-1)
        G_BP.(Slice{i}) = graph(Networkdata_BP.(Slice{i}).opvizmat);
        h=plot(G_BP.(Slice{i}),'NodeLabel',[]);
        highlight(h,Networkdata_BP.(Slice{i}).numbanks+1:(Networkdata_BP.(Slice{i}).numbanks+Networkdata_BP.(Slice{i}).numassets),'NodeColor','r');
        
        %labelnode(h,Networkdata_BP.(Slice{i}).bankarray,Networkdata_BP.(Slice{i}).bankids)
        %labelnode(h,Networkdata_BP.(Slice{i}).numbanks+1:Networkdata_BP.(Slice{i}).numbanks+Networkdata_BP.(Slice{i}).numassets,...
            %Networkdata_BP.(Slice{i}).assetids)
        
% Reorganising to visualise as bipartite graph
        h.XData(1:Networkdata_BP.(Slice{i}).numbanks) = 1;
        %h.YData(1:Networkdata_BP.(Slice{i}).numbanks) = linspace(0.5-(Networkdata_BP.(Slice{i}).numbanks/(2*Networkdata_BP.(Slice{i}).numassets)),...
            %0.5+(Networkdata_BP.(Slice{i}).numbanks/(2*Networkdata_BP.(Slice{i}).numassets)),Networkdata_BP.(Slice{i}).numbanks);
        h.YData(1:Networkdata_BP.(Slice{i}).numbanks) = linspace(0,1,Networkdata_BP.(Slice{i}).numbanks);
        h.XData(Networkdata_BP.(Slice{i}).numbanks+1:(Networkdata_BP.(Slice{i}).numbanks+Networkdata_BP.(Slice{i}).numassets)) = 2;
        h.YData(Networkdata_BP.(Slice{i}).numbanks+1:(Networkdata_BP.(Slice{i}).numbanks+Networkdata_BP.(Slice{i}).numassets))...
            = linspace(0,1,Networkdata_BP.(Slice{i}).numassets);
        title(graphtitle{i},'Interpreter','latex')
        set(gca,'XtickLabel',[],'YtickLabel',[],'XTick',[],'YTick',[]);    
end

set(gcf,'renderer','painters');
hfig = tightfig;
%set(gcf,'Units','Inches');
%pos = get(gcf,'Position');
%set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig,'-dpdf',strcat(fig_output_VN,'overlappingportfolio_evolution.pdf'));

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