function [node] = DANSE_init(node,sim_param,DANSE_param)
% DANSE_INIT - initilization of nodes for DANSE algorithm 
% Syntax:  [node] = DANSE_init(node,sim_param,DANSE_param)
% Inputs:   node            -   structure containing node data
%           DANSE_param     -   DANSE parameters
%           sim_param       -   simulation parameters
%                                                         
% Outputs:  node            -   return structure containing initial cost 
%                               and transmitted signals                        
%
% Other m-files required: filterwithW
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 12-Oct-2015
%------------- BEGIN CODE --------------
% inital cost variable
[node.cost] = deal(0);
%% generate initial broadcast signals
for idx_node = 1:DANSE_param.nb_nodes;    
    for frame_idx = 1:sim_param.ds_idx;
        for ll=1:sim_param.fftL/2+1
            node(idx_node).ZY(:,ll,frame_idx) = ...
                node(idx_node).loc_filt_coeff(:,:,ll)'*node(idx_node).ds_frame(:,ll,frame_idx);
        end
    end
    for frame_idx = 1:sim_param.n_idx;
        for ll=1:sim_param.fftL/2+1
            node(idx_node).ZN(:,ll,frame_idx) = ...
                node(idx_node).loc_filt_coeff(:,:,ll)'*node(idx_node).n_frame(:,ll,frame_idx);
        end
    end
    % broadcast signals in time domain (only used for cost calculation)
    node(idx_node).zy = ...
        filterwithW([node(idx_node).ss_clean],node(idx_node).loc_filt_coeff);
    node(idx_node).zn = ...
        filterwithW([node(idx_node).ss_noise],node(idx_node).loc_filt_coeff);
    
end
%%  calculate initial cost
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
%------------- END OF CODE --------------

