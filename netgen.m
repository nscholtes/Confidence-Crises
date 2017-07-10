function [a,ibn_adjmat] = netgen(n_banks,gamma,a_min,a_max,fileID_IBN,fig_output)

%--------------------------------------------------------------------------

fig_output_IBN   = strcat(fig_output,'Network/');
fig_output_Gephi = strcat(fig_output_IBN,'Gephi/');

% Simulate an undirected network using the 'fitness-based model with mutual
% benefit' - Comprises the following steps:

%%%  - Drawing node fitness distribution from a truncated power law
%%%  - Populating the adjacency matrix
%%%  - Ensuring that the network does not feature isolated nodes
%%%  - Output to GEPHI for network visualisation

%%% Copyright (C) Nicolas K. Scholtes, 2016
%%% Distributed under GPL v3.0

%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
%% Drawing node fitness distribution from a truncated power law
%-------------------------------------------------------------------------

% Monte Carlo simulation with $n_draws$ draws from the distribution using
% the inverse transform sampling method

n_draws = 100000;
n_replicates = n_banks;
n_notconnected = 0;
rnd_LB = 1;
rnd_UB = n_banks/5;

% User-specified minimum degree of all nodes in the network

mindeg = 1;

a_minvec  = ones(n_draws,1)'*a_min;
a_maxvec  = ones(n_draws,1)'*a_max;
a_diffvec = a_maxvec.^(1-gamma) - a_minvec.^(1-gamma);

u = rand(1,n_draws);
all_draws = (a_minvec.^(1-gamma)+ u.*(a_diffvec)).^(1/(1-gamma));

probability_matrix = zeros(n_banks);
adjacency_matrix = zeros(n_banks,n_banks,n_replicates);

% Logarithmic binning to obtain final fitness vector

[a,a_freq] = lnbin(all_draws,n_banks);
a = sort(a,'descend');
a_max = max(a);

figure
hist(a);
title('Empirical distribution of node fitness values')
xlabel('a')
ylabel('f(a)')
set(gcf,'renderer','painters');
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'-dpdf',strcat(fig_output_IBN,'sizedist.pdf'));

%-------------------------------------------------------------------------
%% Populating the adjacency matrix
%-------------------------------------------------------------------------

% Computing probability that two nodes form an undirected edge based on
% their relative fitness

for i=2:n_banks
    j=1;
    while j<i           
        probability_matrix(i,j) = (a(i)/(a_max))*(a(j)/a_max);
        j=j+1;       
    end
end

% Monte Carlo simulation of adjacency matrix (populate below main-diagonal
% for improved computational speed)

for k=1:n_replicates
  rand_mat(:,:,k) = tril(rand(n_banks));
    for i=1:n_banks
        j=1;
        while j<=i
            if rand_mat(i,j,k) < probability_matrix(i,j)                
                adjacency_matrix(i,j,k) = 1;                
            else
                adjacency_matrix(i,j,k) = 0;
            end              
        j=j+1;
        end
    end                 
end    

% Reflect about main-diagonal to obtain final adjacency matrix for the simple network + contingency for isolated nodes (rerun probability
% function for 0-degree nodes until a link is formed and the replicate is connected

for k=1:n_replicates
        adjacency_matrix(:,:,k) = tril(adjacency_matrix(:,:,k),0)+...
        triu(adjacency_matrix(:,:,k)',1);
    
    if chksymmetry(adjacency_matrix(:,:,k))
        fprintf(fileID_IBN,'Adjacency matrix is symmetric. Network is undirected\r\n');
   else
        fprintf(fileID_IBN,'Error! Adjacency matrix is not symmetric\r\n');
    end
    
    [S,count_condviol] = minimum_degree(adjacency_matrix(:,:,k),mindeg);
    
    if  ~S
            n_notconnected = n_notconnected+1;
            fprintf(fileID_IBN,'Replicate matrix %d has %d isolated nodes!\r\n',k,count_condviol);
    elseif S
            fprintf(fileID_IBN,'Replicate matrix %d satisfies the node minimum degree requirement!\r\n',k);
    end
end

fprintf(fileID_IBN,'Percent disconnected replicates in first simulation: %d\r\n', (n_notconnected/n_replicates)*100);

%-------------------------------------------------------------------------
%% Ensuring that the network does not feature isolated nodes
%-------------------------------------------------------------------------

for k=1:n_replicates
        while ~minimum_degree(adjacency_matrix(:,:,k),mindeg)
           
            [~,isolatednode_id_temp] = ismember(adjacency_matrix(:,:,k),zeros(1,n_banks),'rows');
            isolatednode_id = find(isolatednode_id_temp);
            n_isolatednodes = numel(isolatednode_id);
            connect_isolatednode = zeros(1,n_isolatednodes);
                        
            for i =1:n_isolatednodes
                connect_isolatednode(i) = round((rnd_UB-rnd_LB).*rand + rnd_LB);
                adjacency_matrix(isolatednode_id(i),connect_isolatednode(i),k)  = 1;
                adjacency_matrix(connect_isolatednode(i),isolatednode_id(i),k) = 1;
            end
        end
    if minimum_degree(adjacency_matrix(:,:,k),mindeg)
        fprintf(fileID_IBN,'Current disconnected replicate matrix is now connected!');
    else
        fprintf(fileID_IBN,'Still not connected...');
    end
    clearvars isolatenode_id isolatednode_id_temp
end

%-------------------------------------------------------------------------
%% Obtaining  node and undirected edge list (arguments for Gephi network visualization)
%-------------------------------------------------------------------------

% EDGE LIST

for k = 1:n_replicates
    
    [degree,~,~] = degrees(adjacency_matrix(:,:,k));
    totdeg_vec(k) = sum(degree);
end

rand_replicate = round((n_banks-1).*rand + 1);

edge_list =adj2edgeL(adjacency_matrix(:,:,rand_replicate));
header_edges = {'Source' 'Target'};
csvwrite_alt(strcat(fig_output_Gephi,'ibn_edgelist_',datestr(datetime('today')),'.csv'), edge_list, header_edges);

ibn_adjmat = adjacency_matrix(:,:,rand_replicate);

% NODE LIST

nodeid = linspace(1,n_banks,n_banks)';
nodesize = a;
header_nodes = {'Id' 'Size'};
node_data = [nodeid nodesize];
csvwrite_alt(strcat(fig_output_Gephi,'ibn_nodelist_',datestr(datetime('today')),'.csv'), node_data, header_nodes);


end