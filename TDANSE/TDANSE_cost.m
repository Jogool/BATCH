function node = TDANSE_cost(node,sim_param,DANSE_param)
% TDANSE_COST - cost at node for TDANSE algorithm
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
% Oct. 2015; Last revision: 03-Nov-2015
%------------- BEGIN CODE --------------
% inital cost variable
[node.cost] = deal(0);
nb_ds = DANSE_param.desired_sources;
for idx_node = 1:sim_param.nb_nodes
    % adjacency matrix for minimum spanning tree
    A_mst(idx_node,node(idx_node).tree_conn) = 1;
    [rows,cols] = find(A_mst);
    % create sparse matrix
    UG = sparse(rows,cols,1,sim_param.nb_nodes,sim_param.nb_nodes);
end
for idx_node =  1:sim_param.nb_nodes;
    
    idx_nodes = find(idx_node ~= 1:sim_param.nb_nodes);
    for idx = idx_nodes
        node(idx_node).WqqGkq(idx).coeff = node(idx).loc_filt_coeff;
    end
    % construct [Wkk WkkGkq....]
    for idx = idx_nodes
        % find path between current node and node we wish to update
        [~,path,~] = graphshortestpath(UG,idx_node,idx);
        for idx_path = length(path):-1:2
            for ll=1:sim_param.fftL/2+1
                node(idx_node).WqqGkq(idx).coeff(:,:,ll) = ...
                    node(idx_node).WqqGkq(idx).coeff(:,:,ll)*...
                    node(path(idx_path-1)).gkq(path(idx_path)).coeff(:,:,ll);
            end
        end
    end
    % stack filter coefficients
    W = cat(1,cat(1, node(idx_node).loc_filt_coeff),...
        cat(1,node(idx_node).WqqGkq(idx_nodes).coeff));
    idx = [idx_node idx_nodes];
    % estimate desired signal
    dest = filterwithW([node(idx).ss_clean]+[node(idx).ss_noise],W);
    % calculate cost
    node(idx_node).cost = ...
        norm(node(idx_node).ss_clean(:,1:nb_ds) - dest)^2;
end
%------------- END OF CODE --------------
