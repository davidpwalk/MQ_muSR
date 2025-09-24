function [fdy, fdx] = getmultspec(tdy,tdx,Nmult,win)
% function that calculates the spectrum of a 2D time series [fdy, fdx] =
% getmultspec(tdy,tdx,Nmult,win) The input tdy has the samples along its
% first dimension, whereas the second dimension contains multiple data
% series Nmult is length the multiplier for zero padding. win is an
% indicator for window usage. set equal to one for a normalized cheby
% window symmetrically applied (echoe et al), set 2 for asymmetric (FIDs et
% al) fdy is the complex spectrum, normalized by the length fdx is the
% frequency axis

% window time domain data
siz = size(tdy);
if win == 1
    win = chebwin(siz(1));
%     winM = repmat(win,1,siz(2));
    tdw = bsxfun(@times,tdy,win/mean(win));
elseif win == 2
    win_2 = ifftshift(chebwin(siz(1)*2));
	win = win_2(1:siz(1));
    tdw = bsxfun(@times,tdy,win/mean(win));
else
    tdw = tdy;
end

%this part of the code generates that frequency axis
N=siz(1)*Nmult;
if mod(N,2)==0
    k=-N/2:N/2-1; % N even
else
    k=-(N-1)/2:(N-1)/2; % N odd
end
T=N*(tdx(2)-tdx(1));
fdx=k/T;  %the frequency axis
fdx = transpose(fdx);


%takes the fft of the signal
fdy=fft(tdw,N)/siz(1); % normalize the data
fdy=fftshift(fdy,1); %shifts the fft data such that it is centered