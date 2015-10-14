function [node] = centralized_batch(node,DANSE_param,sim_param)
% CENTRALIZED_BATCH - perform centralized batch mode simulation 
% Syntax:  [node] = centralized_batch(node,DANSE_param,sim_param)
% Inputs:   node            -   structure containing node data
%           DANSE_param     -   DANSE parameters
%           sim_param       -   simulation parameters
%                                                         
% Outputs:  node            - structure containing node data, returned with
%                               centralized cost and centralized filters
%
% Other m-files required: filterwithW
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 11-Oct-2015
%------------- BEGIN CODE --------------
% initialize correlation matrices
Ryy = zeros(DANSE_param.nb_nodes*DANSE_param.sensors,...
    DANSE_param.nb_nodes*DANSE_param.sensors,...
    sim_param.fftL/2+1);
Rnn = Ryy;

% collect desired source frames from all nodes
ds_frames = cat(1,node.ds_frame);
% collect noise frames from all nodes
n_frames =  cat(1,node.n_frame);

for frame_idx = 1:sim_param.ds_idx;
    for ll=1:sim_param.fftL/2+1
        Y = ds_frames(:,ll,frame_idx);
        Ryy(:,:,ll)=Ryy(:,:,ll)+Y*Y';
    end
end

for frame_idx = 1:sim_param.n_idx;
    for ll=1:sim_param.fftL/2+1
        Y = n_frames(:,ll,frame_idx);
        Rnn(:,:,ll)=Rnn(:,:,ll)+Y*Y';
    end
end

Ryy = Ryy /  sim_param.ds_idx;
Rnn = Rnn / sim_param.n_idx;
%!!! change the following code to estimate Rxx
Rxx = Ryy;
%% filter and cost calculation
% filter
for ll=1:sim_param.fftL/2+1
    W(:,:,ll) = (Rnn(:,:,ll)+Rxx(:,:,ll)) \ Rxx(:,:,ll);
end
% calculate cost and store centralized filter
% idx_sen : sensor index
idx_sen = 1;
% idx_node : node index
for idx_node = 1:DANSE_param.nb_nodes

    dest = ...
        filterwithW([node.ss_clean] + [node.ss_noise],...
        W(:,idx_sen:idx_sen+DANSE_param.desired_sources-1,:));
    
    node(idx_node).cent_cost = ...
        norm(node(idx_node).ss_clean(:,1:DANSE_param.desired_sources) - dest)^2;
    node(idx_node).cent_filt = ...
       W(:,idx_sen:idx_sen+DANSE_param.desired_sources-1,:);
    
    idx_sen = idx_sen+node(idx_node).sensors;
end
%------------- END OF CODE --------------