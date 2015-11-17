function [node] = TDANSE_ff(node,node_update,sim_param,DANSE_param)
% TDANSE_FF - given a sink (root) node and pre-existing tree, find the
%               data flow toward the sink (root) node
% Syntax:  node = TDANSE_ff(node,sim_param,DANSE_param)
% Inputs:   node            -   structure containing node data
%           node_update     -   updating node (root) node of the tree
%           sim_param       -   simulation parameters
%           DANSE_param     -   DANSE parameters
%
% Outputs:  node            -   return structure containing FF signals
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 03-Nov-2015
%------------- BEGIN CODE --------------
nb_ds = DANSE_param.desired_sources;
[node.ff_update] = deal(0);   % flag if node has transmitted its ff signal

% leaf node FF signals are equivalent to DANSE signals
for idx_node = find(cellfun(@(x) numel(x), {node.tree_conn}) == 1)
    if ~eq(idx_node,node_update)

        node(idx_node).ff_ZY = node(idx_node).loc_ZY;
        node(idx_node).ff_ZN = node(idx_node).loc_ZN;
        
        % set FF flag
        node(idx_node).ff_update = 1;
    end
end

% find all nodes who have performed a FF update
node_ff_update = numel(find(cat(1,node.ff_update)));

while lt(node_ff_update,DANSE_param.nb_nodes-1) % can skip the root node, hence - 1
    % find all nodes who have not performed a fusion flow update
    idx_ff = find(~cat(1,node.ff_update));
    % remove the root node as it will never generate a ff signal
    idx_ff(idx_ff == node_update) = [];
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

            % gather all updated neighbor signals
            ZY = cat(1,node(nbrs_updated).ff_ZY);
            ZN = cat(1,node(nbrs_updated).ff_ZN);
            
            % gather all Gkq coefficients to apply to neighbor signals
            Gkq_coeff = cat(1,node(idx_node).gkq(nbrs_updated).coeff);
            
            Gkq_size = size(Gkq_coeff,1);
            if eq(Gkq_size,1)
                for idx_z = 1:nb_ds
                    for ll=1:sim_param.fftL/2+1
                        node(idx_node).ff_ZY(idx_z,ll,:) = ...
                            node(idx_node).loc_ZY(idx_z,ll,:)+ ...
                            bsxfun(@times,conj(Gkq_coeff(:,idx_z,ll)),ZY(:,ll,:));
                        
                         node(idx_node).ff_ZN(idx_z,ll,:) = ...
                            node(idx_node).loc_ZN(idx_z,ll,:)+ ...
                            bsxfun(@times,conj(Gkq_coeff(:,idx_z,ll)),ZN(:,ll,:));
                    end
                end
            else
                for idx_z = 1:nb_ds
                    for ll=1:sim_param.fftL/2+1
                        node(idx_node).ff_ZY(idx_z,ll,:) = ...
                            node(idx_node).loc_ZY(idx_z,ll,:)+ ...
                            sum(bsxfun(@times,conj(Gkq_coeff(:,idx_z,ll)),ZY(:,ll,:)));
                        
                        node(idx_node).ff_ZN(idx_z,ll,:) = ...
                            node(idx_node).loc_ZN(idx_z,ll,:)+ ...
                            sum(bsxfun(@times,conj(Gkq_coeff(:,idx_z,ll)),ZN(:,ll,:)));
                    end
                end
            end
 
            % set ff flag
            node(idx_node).ff_update = 1;
        end
    end
    % find all nodes who have performed a fusion flow update
    node_ff_update = numel(find(cat(1,node.ff_update)));
end
%------------- END OF CODE --------------