function node = TIDANSE_cost(node,sim_param,DANSE_param)
% TIDANSE_COST - cost at node for TIDANSE algorithm
% Syntax:  [node] = DANSE_init(node,sim_param,DANSE_param)
% Inputs:   node            -   structure containing node data
%           DANSE_param     -   DANSE parameters
%           sim_param       -   simulation parameters
%                                                         
% Outputs:  node            -   initial cost and transmitted signals                        
%
% Other m-files required: filterwithW
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 02-Nov-2015
%------------- BEGIN CODE --------------
% inital cost variable
[node.cost] = deal(0);
nb_ds = DANSE_param.desired_sources;
%%  calculate initial cost
for idx_node =  1:sim_param.nb_nodes;
    
    idx_nodes = find(idx_node ~= 1:sim_param.nb_nodes);
    % construct [Wkk WkkGkq....]
    for idx = idx_nodes
        for ll=1:sim_param.fftL/2+1
            node(idx_node).WqqGkq(idx).coeff(:,:,ll) = ...
                node(idx).P(:,:,ll)...
                *node(idx_node).gkq(1).coeff(:,:,ll);
        end
    end
    
    W = cat(1,cat(1, node(idx_node).loc_filt_coeff),...
        cat(1,node(idx_node).WqqGkq(idx_nodes).coeff));
    % estimated desired signal
    idx = [idx_node idx_nodes];
    dest = filterwithW([node(idx).ss_clean]+[node(idx).ss_noise],W);
    % calculate cost
    node(idx_node).cost = ...
        norm(node(idx_node).ss_clean(:,1:nb_ds) - dest)^2;
end
%------------- END OF CODE --------------
