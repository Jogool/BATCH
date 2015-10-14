function [node] = TDANSE_batch(node,sim_param,DANSE_param,node_update)
% TDANSE_BATCH - Wrapper to run a single iteration of the T-DANSE algorithm 
% Syntax:  [node] = TDANSE_batch(node,sim_param,DANSE_param,node_update)
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
% inital cost variable
[node.cost] = deal(0);
% fusion flow
node = TDANSE_ff(node,node_update,sim_param,DANSE_param);
% diffusion flow
node = TDANSE_df(node,node_update,sim_param);
% filter update
node = TDANSE_filt_update(node,sim_param,DANSE_param,node_update);

% cost calculation
for idx_node = 1:DANSE_param.nb_nodes;  
    idx = node(idx_node).ff_rec;
    % gather local signals and z-signals of other nodes
    y = [node(idx_node).ss_clean + node(idx_node).ss_noise];
    zy = [node(node(idx_node).ff_trans).df(idx_node).zy node(idx).ff_zy];
    zn = [node(node(idx_node).ff_trans).df(idx_node).zn node(idx).ff_zn];
    
    % gather G coefficiens
    Gkq_coeff = ...
        cat(1,[node(idx_node).gkq(node(idx_node).ff_trans).coeff],cat(1,node(idx_node).gkq(idx).coeff));
    
    W = cat(1, node(idx_node).loc_filt_coeff, Gkq_coeff);
    % estimated desired signal
    dest = filterwithW([y zy+zn],W);
    % calculate cost
    node(idx_node).cost = ...
        norm(node(idx_node).ss_clean(:,1:DANSE_param.desired_sources) - dest)^2;
end
%% Update local broadcast signals of update node in time domain (only used for cost calculation)
node(node_update).loc_zy = ...
    filterwithW([node(node_update).ss_clean],node(node_update).loc_filt_coeff);
node(node_update).loc_zn = ...
    filterwithW([node(node_update).ss_noise],node(node_update).loc_filt_coeff);
%------------- END OF CODE --------------
