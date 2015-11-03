function node = TDANSE_ff(node,node_update,sim_param,DANSE_param)
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
[node.ff_update] = deal(0);   % flag if node has transmitted its FF signal
nb_ds = DANSE_param.desired_sources;
for idx_node = 1:sim_param.nb_nodes
    % adjacency matrix for current tree
    A_mst(idx_node,node(idx_node).tree_conn) = 1;
    [rows,cols] = find(A_mst);
    % create sparse matrix
    UG = sparse(rows,cols,1,sim_param.nb_nodes,sim_param.nb_nodes);
end

% find all leaf nodes
idx_nodes = find(cellfun(@(x) numel(x), {node.tree_conn}) == 1);
% remove updating node if a leaf node
idx_nodes(idx_nodes == node_update) = [];

for idx_node = idx_nodes
    % find path from leaf node to root node
    [~,path,~] = graphshortestpath(UG,node_update,idx_node);
    % remove root node
    path(1) = [];
    % leaf node FF signals are equivalent to DANSE signals
    node(idx_node).ff_ZY = node(idx_node).loc_ZY;
    node(idx_node).ff_ZN = node(idx_node).loc_ZN;
    % update FF flag for leaf node
    node(idx_node).ff_update = 1;
    % if the path goes through another node, have to combine FF signals
    if ~eq(length(path),1)
        % combine FF signals from the leaf node to the root node
        for idx_path = length(path):-1:2
            % gather all Gkq coefficients to apply to neighbor signals
            Gkq_coeff = cat(1,node(path(idx_path-1)).gkq(path(idx_path)).coeff);
            % if the node has already updated its FF signal, then do not add
            % the local contribution
            if node(path(idx_path-1)).ff_update
                for idx_z = 1:nb_ds
                    for ll=1:sim_param.fftL/2+1
                        node(path(idx_path-1)).ff_ZY(idx_z,ll,:) = ...
                            node(path(idx_path-1)).ff_ZY(idx_z,ll,:)+ ...
                            sum(bsxfun(@times,conj(Gkq_coeff(:,idx_z,ll)),...
                            node(path(idx_path)).ff_ZY(:,ll,:)));
                        
                        node(path(idx_path-1)).ff_ZN(idx_z,ll,:) = ...
                            node(path(idx_path-1)).ff_ZN(idx_z,ll,:)+ ...
                            sum(bsxfun(@times,conj(Gkq_coeff(:,idx_z,ll)),...
                            node(path(idx_path)).ff_ZN(:,ll,:)));
                    end
                end
            else
                for idx_z = 1:nb_ds
                    for ll=1:sim_param.fftL/2+1
                        node(path(idx_path-1)).ff_ZY(idx_z,ll,:) = ...
                            node(path(idx_path-1)).loc_ZY(idx_z,ll,:)+ ...
                            sum(bsxfun(@times,conj(Gkq_coeff(:,idx_z,ll)),...
                            node(path(idx_path)).ff_ZY(:,ll,:)));
                        
                        node(path(idx_path-1)).ff_ZN(idx_z,ll,:) = ...
                            node(path(idx_path-1)).loc_ZN(idx_z,ll,:)+ ...
                            sum(bsxfun(@times,conj(Gkq_coeff(:,idx_z,ll)),...
                            node(path(idx_path)).ff_ZN(:,ll,:)));
                    end
                end
                % update FF flag
                node(path(idx_path-1)).ff_update = 1;
            end
        end
    end
end
%------------- END OF CODE --------------