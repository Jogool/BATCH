function [node] = DANSE_batch(node,sim_param,DANSE_param,node_update)
% DANSE_BATCH - Wrapper to run a single iteration of the DANSE algorithm 
% Syntax:  [node] = DANSE_batch(node,sim_param,DANSE_param,node_update)
% Inputs:   node            -   structure containing node data
%           DANSE_param     -   DANSE parameters
%           sim_param       -   simulation parameters
%           node_update     -   node to update
%                                                         
% Outputs:  node            -   return structure containing cost, updated
%                                   filter and broadcast signals
%
% Other m-files required: DANSE_filt_update, filterwithW
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 12-Oct-2015
%------------- BEGIN CODE --------------
%% Filter update
% update correlation matrices and filter (only needs to happen on the updating node)
node = DANSE_filt_update(node,sim_param,DANSE_param,node_update);

%% Cost calculation
for idx_node = 1:DANSE_param.nb_nodes;  
    % find all nodes ~= to idx_node
    idx = find(idx_node ~= 1:DANSE_param.nb_nodes);
    % gather all W and G coefficients
    W = cat(1, node(idx_node).loc_filt_coeff, node(idx_node).gkq.coeff);
    % total received signal at node
    y = [node(idx_node).ss_clean + node(idx_node).ss_noise [node(idx).zy]+[node(idx).zn]];
    % estimated desired signal
    dest = filterwithW(y,W);
    % calculate cost
    node(idx_node).cost = ...
        norm(node(idx_node).ss_clean(:,1:DANSE_param.desired_sources) - dest)^2;
end
%% Update broadcast signals of update node in time domain (only used for cost calculation)
    node(node_update).zy = ...
        filterwithW([node(node_update).ss_clean],node(node_update).loc_filt_coeff);
    node(node_update).zn = ...
        filterwithW([node(node_update).ss_noise],node(node_update).loc_filt_coeff);
%------------- END OF CODE --------------
