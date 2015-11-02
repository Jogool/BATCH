function [node] = TIDANSE_batch(node,sim_param,DANSE_param,node_update)
% TIDANSE_BATCH - Wrapper to run a single iteration of the TI-DANSE algorithm 
% Syntax:  [node] = TIDANSE_batch(node,sim_param,DANSE_param,node_update)
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
% Oct. 2015; Last revision: 13-Oct-2015
% !!!!!! the filters converge to the correct values, however, there is a
% problem with using the fftfilt command to calculate the cost.  This
% probably has something to do with applying the filter at a node in the
% frequency domain, converting to time domain, and repeating this process.
% I would look into probably multiplying all of the filters in frequency
% domain somehow and then just apply to these to the local signals at the
% nodes.
%------------- BEGIN CODE --------------
% inital cost variable
[node.cost] = deal(0);

% filter update
node = TIDANSE_filt_update(node,sim_param,DANSE_param,node_update);
% find cost
node = TIDANSE_cost(node,sim_param,DANSE_param);

% %% Cost calculation
% zy = [node.loc_zy];
% zn = [node.loc_zn];
% 
% for idx_DANSE = 1:DANSE_param.desired_sources;
%     zy_sum(:,idx_DANSE) = ...
%         sum(zy(:,idx_DANSE:DANSE_param.desired_sources:end),2);
%     zn_sum(:,idx_DANSE) = ...
%         sum(zn(:,idx_DANSE:DANSE_param.desired_sources:end),2);
% end
% 
% for idx_node = 1:DANSE_param.nb_nodes;  
%     zy_loc = zy_sum - node(idx_node).loc_zy;
%     zn_loc = zn_sum - node(idx_node).loc_zn;
% 
% 
%     W = cat(1, node(idx_node).loc_filt_coeff, node(idx_node).gkq.coeff);
%     y = [node(idx_node).ss_clean + node(idx_node).ss_noise zy_loc+zn_loc];
% 
%     
%     % estimated desired signal
%     dest = filterwithW(y,W);
%     % calculate cost
%     node(idx_node).cost = ...
%         norm(node(idx_node).ss_clean(:,1:DANSE_param.desired_sources) - dest)^2;
% end
% %% Update local broadcast signals of update node in time domain (only used for cost calculation)
% node(node_update).loc_zy = ...
%    filterwithW(node(node_update).ss_clean,node(node_update).P);
% node(node_update).loc_zn = ...
%    filterwithW(node(node_update).ss_noise,node(node_update).P);
%------------- END OF CODE --------------
