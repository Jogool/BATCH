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
% Other m-files required: TIDANSE_filt_update, TIDANSE_cost
% Subfunctions: none
% MAT-files required: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 02-Nov-2015
%------------- BEGIN CODE --------------
% filter update
node = TIDANSE_filt_update(node,sim_param,DANSE_param,node_update);
% find cost
node = TIDANSE_cost(node,sim_param,DANSE_param);
%------------- END OF CODE --------------
