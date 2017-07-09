function [S,count_condviol] = minimum_degree(adj,mindeg)
 
% Check to see that all edges have a minimum user-specified degree
 
nodes_condsat  = find(sum(adj) >= mindeg);
nodes_condviol = find(sum(adj) < mindeg);
count_condviol = numel(nodes_condviol);
 
% S is a binary parameter that is true only if there are NO isolated nodes
% and false if there exists at least one
 
if count_condviol ~= 0   
    S= false;   
else
    S = true;
end
 
end