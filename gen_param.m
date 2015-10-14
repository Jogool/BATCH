function [sim_param,DANSE_param] = gen_param(source_sig,noise_sig)
% GEN_PARAM - generate simulation and DANSE parameters
% Syntax:  [] = frame_signal(node,sim_param)
% generates frames of data for batch algorithms
% Inputs:   source_sig      -   single channel desired source signal
%           noise_sig       -   additive noise (can be multichannel)
%                                                         
% Outputs:  DANSE_param     -   structure with DANSE parameters
%           sim_param       -   structure with simulation parameters
%
% Other m-files: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 13-Oct-2015
%------------- BEGIN CODE --------------
%% simulation parameters
sim_param.L = [5 4 6];                                                      % room dimensions
sim_param.fftL = 128;                                                        % FFT block size
sim_param.length_of_signal = length(source_sig);                            % length of signal
sim_param.plot = 1;                                                         % 1(0) - show (do not show) network
sim_param.nb_nodes = 5;
sim_param.ds_idx = 0;                                                       % counter for desired signal activity
sim_param.n_idx = 0;                                                        % counter for interferer activity
sim_param.overlap = floor(sim_param.fftL/2);                                % frame overlap
sim_param.nb_frames = ...
    floor((sim_param.length_of_signal-sim_param.fftL)/ ...
    (sim_param.fftL-sim_param.overlap) + 1);                                 % number of frames

%% DANSE parameters
DANSE_param.desired_sources = size(source_sig,2);                           % number of desired sources (also dimension of DANSE)
DANSE_param.noise_sources = size(noise_sig,2);                              % number of correlated noise sources            
DANSE_param.nb_nodes = sim_param.nb_nodes;                                  % number of nodes
DANSE_param.sensors = DANSE_param.desired_sources + 1;                      % number of sensors per node (assumed same across all nodes compression 
DANSE_param.DANSE_update = 0;   
DANSE_param.thresh = 1e-5;                                                  % stopping threshold for DANSE
DANSE_param.max_iter = 250;                                                  % max iterations for DANSE
%------------- END OF CODE --------------
