function [node] = TDANSE_df(node,idx_node,sim_param)
% TDANSE_DF - given a sink (root) node and pre-existing tree, find the
%                   data flow away from the root node (recursively)
%
% Syntax:  [node] = TDANSE_df(node,sim_param,DANSE_param)
% Inputs:   node            -   structure containing node data
%           root            -   root (sink) node of the tree
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
for idx_ffrec = node(idx_node).ff_rec
    idx = node(idx_node).ff_rec;
    idx = idx(idx~= idx_ffrec);
    
    if isempty(idx) % single neighbor in the diffusion flow
        
        % Note that there are times when the root node does not have any
        % other signals besides its own to broadcast in the diffusion flow
        %, i.e., a line topology when the root is at one of the ends.
        % Therefore this if statement catches this condition
        if isempty(node(idx_node).ff_trans)
            % update df in frequency domain
            node(idx_node).df(idx_ffrec).ZY = node(idx_node).loc_ZY;
            node(idx_node).df(idx_ffrec).ZN = node(idx_node).loc_ZN;
            
            % update df in time domain
            node(idx_node).df(idx_ffrec).zy =  node(idx_node).loc_zy;
            node(idx_node).df(idx_ffrec).zn =  node(idx_node).loc_zn;
            
        else
            ZY = node(node(idx_node).ff_trans).df(idx_node).ZY;
            ZN = node(node(idx_node).ff_trans).df(idx_node).ZN;
            
            Gkq_coeff = node(idx_node).gkq(node(idx_node).ff_trans).coeff;
            
            % update df in frequency domain
            for frame_idx = 1:sim_param.ds_idx;
                for ll=1:sim_param.fftL/2+1
                    node(idx_node).df(idx_ffrec).ZY(:,ll,frame_idx) = ...
                        node(idx_node).loc_ZY(:,ll,frame_idx) + ...
                        Gkq_coeff(:,:,ll)'*ZY(:,ll,frame_idx);
                end
            end
            for frame_idx = 1:sim_param.n_idx;
                for ll=1:sim_param.fftL/2+1
                    node(idx_node).df(idx_ffrec).ZN(:,ll,frame_idx) = ...
                        node(idx_node).loc_ZN(:,ll,frame_idx) + ...
                        Gkq_coeff(:,:,ll)'*ZN(:,ll,frame_idx);
                end
            end
            
            % update df in time domain
            node(idx_node).df(idx_ffrec).zy = node(idx_node).loc_zy + ...
                filterwithW([node(node(idx_node).ff_trans).df(idx_node).zy],Gkq_coeff);
            node(idx_node).df(idx_ffrec).zn = node(idx_node).loc_zn + ...
                filterwithW([node(node(idx_node).ff_trans).df(idx_node).zn],Gkq_coeff);
            
        end
    else % multiple neighbors in the diffusion flow
        if isempty(node(idx_node).ff_trans) % root node condition
            ZY = cat(1,node(idx).ff_ZY);
            ZN = cat(1,node(idx).ff_ZN);
            
            zy = [node(idx).ff_zy];
            zn = [node(idx).ff_zn];
            
        else
            ZY = cat(1,node(node(idx_node).ff_trans).df(idx_node).ZY,node(idx).ff_ZY);
            ZN = cat(1,node(node(idx_node).ff_trans).df(idx_node).ZN,node(idx).ff_ZN);
            
            zy = [node(node(idx_node).ff_trans).df(idx_node).zy node(idx).ff_zy];
            zn = [node(node(idx_node).ff_trans).df(idx_node).zn node(idx).ff_zn];
        end
        
        % get gkq coefficeints used for fusion
        Gkq_coeff = ...
            cat(1,node(idx_node).gkq(node(idx_node).ff_trans).coeff, cat(1,node(idx_node).gkq(idx).coeff));
        
        % update df in frequency domain
        for frame_idx = 1:sim_param.ds_idx;
            for ll=1:sim_param.fftL/2+1
                node(idx_node).df(idx_ffrec).ZY(:,ll,frame_idx) = ...
                    node(idx_node).loc_ZY(:,ll,frame_idx) + ...
                    Gkq_coeff(:,:,ll)'*ZY(:,ll,frame_idx);
            end
        end
        for frame_idx = 1:sim_param.n_idx;
            for ll=1:sim_param.fftL/2+1
                node(idx_node).df(idx_ffrec).ZN(:,ll,frame_idx) = ...
                    node(idx_node).loc_ZN(:,ll,frame_idx) + ...
                    Gkq_coeff(:,:,ll)'*ZN(:,ll,frame_idx);
            end
        end
        
        % update df in time domain
        node(idx_node).df(idx_ffrec).zy = node(idx_node).loc_zy + ...
            filterwithW(zy,Gkq_coeff);
        node(idx_node).df(idx_ffrec).zn = node(idx_node).loc_zn + ...
            filterwithW(zn,Gkq_coeff);
    end
    
    % recursive call, go all the way to the leaf nodes of the branch
    node = TDANSE_df(node,idx_ffrec,sim_param); 
end
%------------- END OF CODE --------------