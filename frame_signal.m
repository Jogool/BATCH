function [node,sim_param] = frame_signal(node,sim_param)
% FRAME_SIGNAL - convert time domain signal to frequency domain blocks
% Syntax:  [] = frame_signal(node,sim_param)
% generates frames of data for batch algorithms
% Inputs:   node            -  generated nodes in a structure
%           sim_param       -   structure with simulation parameters
%                                                         
% Outputs:  frames          - structure with desired source frames
%                               (frames.ds) and noise frames (frames.n) and
%                               VAD decision (frames.VAD)
%           sim_param       - returns sim param with count of desired
%                               signal and noise frames
%
% Other m-files: none
%
% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Oct. 2015; Last revision: 11-Oct-2015
%------------- BEGIN CODE --------------
for frame_idx = 1:sim_param.nb_frames
    fs = (frame_idx-1)*(sim_param.fftL-sim_param.overlap ) + 1;     % start of frame
    fe = fs+sim_param.fftL-1;                                       % end of frame
    
    VAD_dec = sumsqr(node(1).ss_clean(fs:fe,:)) > .5;
    
    if VAD_dec
        % sim_param.ds_idx: counter for number of desired signal frames
        sim_param.ds_idx = sim_param.ds_idx + 1;
        for node_idx = 1:size(node,2)
            Y = fft(node(node_idx).ss_clean(fs:fe,:)).';
            Y = Y(:,1:sim_param.fftL/2+1);
            node(node_idx).ds_frame(:,:,sim_param.ds_idx) = Y;
        end
    else
        % sim_param.n_idx: counter for number of noise signal frames
        sim_param.n_idx = sim_param.n_idx + 1;
        for node_idx = 1:size(node,2)
            Y = fft(node(node_idx).ss_noise(fs:fe,:)).';
            Y = Y(:,1:sim_param.fftL/2+1);
            node(node_idx).n_frame(:,:,sim_param.n_idx) = Y;
        end
    end
end
%------------- END OF CODE --------------