function [node] = DANSE_batch(node,sim_param,DANSE_param,node_update)
% DANSE_BATCH - Wrapper to run a single iteration of the DANSE algorithm 
% Syntax:  [node] = DANSE_batch(node,sim_param,DANSE_param,node_update)
% Inputs:   node            -   structure containing node data
%           DANSE_param     -   DANSE parameters
%           sim_param       -   simulation parameters
%           node_update     -   node to update
%                                                         
% Outputs:  node            -   return structure containing cost, updated
%                                   filter and broadcast signals
%
% Other m-files required: DANSE_filt_update, DANSE_cost
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 01-Nov-2015
%------------- BEGIN CODE --------------
%% Filter update
% update correlation matrices and filter (only needs to happen on the updating node)
node = DANSE_filt_update(node,sim_param,DANSE_param,node_update);
node = DANSE_cost(node,sim_param,DANSE_param); 
%------------- END OF CODE --------------
