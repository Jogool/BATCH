function [z,w] = filterwithW(y,W)
%NOTICE: this function also compensates for delay introduces by filtering.
%INPUTS:
%   y     = the signals to be filtered 
%   W     = the filters in the frequency domain, 
%            rows : mic, columns : signal subspace dimension, plate : frequency
% Q : desired dimeinons
% L : number of frequency bnds
% M : number of microphones

L = size(W,3)-1;
Q = size(W,2);
M = size(W,1);
for l=2:L
% Conjugate symmetry
    W(:,:,2*L-l+2) = conj(W(:,:,l));
end
% Calculate output signal
z = zeros(size(y,1),Q);
for q=1:Q
    for m = 1:M
    w(:,m) = ifftshift(real(ifft(conj(W(m,q,:)))));  % conj because of W^H, 
    z(:,q) = z(:,q) + fftfilt(w(:,m),y(:,m));
    end
end

%delay compensation
z=[z(L+1:end,:);zeros(L,Q)];
