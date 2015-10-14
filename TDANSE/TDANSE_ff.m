function [node] = TDANSE_ff(node,root,sim_param,DANSE_param)
% TDANSE_FF - given a sink (root) node and pre-existing tree, find the 
%               data flow toward the sink (root) node
% 
% Syntax:  [node] = TDANSE_ff(node,sim_param,DANSE_param)
% Inputs:   node            -   structure containing node data
%           root            -   root (sink) node of the tree
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
[node.ff_trans] = deal([]);   % node k transmits to this node during the ff (should always be a single node)
[node.ff_rec] = deal([]);     % node k receives these signals during the ff
[node.ff_update] = deal(0);   % flag if node has transmitted its ff signal

% if node only has one connection and is not the root node, then it can 
% immediately transmit its ff signal
for idx_node = find(cellfun(@(x) numel(x), {node.tree_conn}) == 1)
    if ~eq(idx_node,root)
        % which node the current node transmits to during the ff
        node(idx_node).ff_trans = node(idx_node).tree_conn;
        % the node that receives the ff from the node
        node(node(idx_node).tree_conn).ff_rec = ...
            sort([node(node(idx_node).tree_conn).ff_rec idx_node]);
        
        % update ff in frequency domain 
        node(idx_node).ff_ZY = node(idx_node).loc_ZY;
        node(idx_node).ff_ZN = node(idx_node).loc_ZN;

        % update ff in time domain 
        node(idx_node).ff_zy =  node(idx_node).loc_zy;
        node(idx_node).ff_zn =  node(idx_node).loc_zn;
        
        % set ff flag
        node(idx_node).ff_update = 1; 
    end
end

% find all nodes who have performed a fusion flow update
node_ff_update = numel(find(cat(1,node.ff_update)));   

while lt(node_ff_update,DANSE_param.nb_nodes-1) % can skip the root node, hence - 1
    % find all nodes who have not performed a fusion flow update
    idx_ff = find(~cat(1,node.ff_update));  
    % remove the root node as it will never generate a ff signal
    idx_ff(idx_ff == root) = [];
    for idx_node = idx_ff'

        % find neighbors who have transmitted a ff signal
        idx = find([node(node(idx_node).tree_conn).ff_update]);
        nbrs_updated = node(idx_node).tree_conn(idx);
        nb_nbrs_updated = numel(nbrs_updated);
        
        % find neighbors who have not transmitted a ff signal
        idx = find(~[node(node(idx_node).tree_conn).ff_update]);
        non_update_neighbors = node(idx_node).tree_conn(idx);
        
        % if all neighbors except 1 have performed a fusion flow
        % update, the node can generate its fusion flow signal
        if eq(nb_nbrs_updated,numel(node(idx_node).tree_conn)-1)    
            
            % update node list of node that will recieve ff signal from
            % node_idx
            idx = sort([node(non_update_neighbors).ff_rec idx_node]);
            node(non_update_neighbors).ff_rec = idx;
            
            % node that node_idx is transmitting to
            node(idx_node).ff_trans = non_update_neighbors;
            
            % gather all updated neighbor signals
            ZY = cat(1,node(nbrs_updated).ff_ZY);
            ZN = cat(1,node(nbrs_updated).ff_ZN);
            
            % gather all Gkq coefficients to apply to neighbor signals
            Gkq_coeff = cat(1,node(idx_node).gkq(nbrs_updated).coeff);
            
            % update ff in frequency domain
            for frame_idx = 1:sim_param.ds_idx;
                for ll=1:sim_param.fftL/2+1
                    node(idx_node).ff_ZY(:,ll,frame_idx) = ...
                        node(idx_node).loc_ZY(:,ll,frame_idx) + ...
                        Gkq_coeff(:,:,ll)'*ZY(:,ll,frame_idx);
                end
            end
            for frame_idx = 1:sim_param.n_idx;
                for ll=1:sim_param.fftL/2+1
                    node(idx_node).ff_ZN(:,ll,frame_idx) = ...
                        node(idx_node).loc_ZN(:,ll,frame_idx) + ...
                        Gkq_coeff(:,:,ll)'*ZN(:,ll,frame_idx);
                end
            end
            
            % update ff in time domain
            node(idx_node).ff_zy = node(idx_node).loc_zy + ...
                filterwithW([node(nbrs_updated).ff_zy],Gkq_coeff);
            node(idx_node).ff_zn = node(idx_node).loc_zn + ...
                filterwithW([node(nbrs_updated).ff_zn],Gkq_coeff);
            
            % set ff flag
            node(idx_node).ff_update = 1;
        end
    end
    % find all nodes who have performed a fusion flow update
    node_ff_update = numel(find(cat(1,node.ff_update)));
end
%------------- END OF CODE --------------