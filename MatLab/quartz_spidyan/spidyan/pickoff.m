function [t] = pickoff(beta, dfrq, nu1, index)
% [t] = balk(beta,dfrq,nu1,index)
% 
% calculates correct pulse duration for a given nu1 and bandwidth
% 
% for further details see
% A. Schweiger, G. Jeschke: Principles of pulse electron paramagnetic
% resonance, Oxford University Press, 2001, p. 137
% 
% Input:
% beta    flip angle given in rad
% dfrq       bandwidth in GHz
% nu1     wave amplitude [MHz] 
% index   number of event pulse experiment
% 
% Output:
% t     calculated pulselength [ns]
%
% balk checks if the conditions for fast passage are fullfilled and will
% print a warning if not.
% 
% S. Pribitzer, 2014

omega1=nu1*2*1e6*pi; % converts irradiation strength into rad
omega=abs(dfrq)*2*1e9*pi; % converts sweeprange to rad
alpha=pi^3*omega1^2/(2*3.01*beta^2); % calculates sweep rate

t=omega/alpha*1e9; % computes required pulse length

if alpha/omega1<0.1
    fprintf('Event %1.0f: Fast passage condition is not satisfied.\n', index);
elseif alpha/omega>0.1 && alpha/omega<1
    fprintf('Event %1.0f: Fast passage condition may be satisfied.\n', index);
else
    fprintf('Event %1.0f: Fast passage condition satisfied.\n', index);
end

end

