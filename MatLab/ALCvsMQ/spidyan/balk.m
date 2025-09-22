function [nu1] = balk(beta, dfrq, t, index)
% [nu1]=balk(beta,dfrq,t,index)
% 
% calculates wave amplitude for a pulse from flip angle 
% duration and bandwidth with Landau Zener Formula
% 
% for further details see
% 1.Jeschke, G., Pribitzer, S. & Doll, A. Coherence Transfer by Passage 
% Pulses in Electron Paramagnetic Resonance Spectroscopy. 
% J. Phys. Chem. B 119, 13570–13582 (2015).
% 
% Input:
% beta    flip angle in rad
% dfrq    bandwidth in GHz
% t       pulselength in ns
% index   number of event pulse experiment
%
% Output:
% nu1     calculated wave amplitude in MHz
%
% balk checks if the conditions for fast passage are fulfilled and will
% print a warning if not.
% 
% S. Pribitzer, 2015


t=t*1e-9; % converts nanoseconds to seconds

if dfrq~= 0
    omega=abs(dfrq)*2*1e9*pi; % sweeprange in rad
       
    alpha=omega/t; % sweeprate in rad s^-1
    
    Q = 2*log(2/(cos(beta)+1))/pi;
    
    if Q > 10 % if Q is larger than 10 (inf in case beta = pi), it is set 
              % to a reasonable high value of 5, which provides adiabatic
              % passage.
        Q = 5;
    end
        
    nu1=sqrt(Q*alpha)/2/pi/1e6;
    
   
    if Q<1
        fprintf('Event %1.0f: Fast passage condition may be satisfied, Q = %1.2f.\n' , index,Q);
    elseif Q>5 && alpha/omega<1
        fprintf('Event %1.0f: Adiabatic passage, Q = %1.2f.\n' , index,Q);
    else
        fprintf('Event %1.0f: Fast passage, Q = %1.2f.\n' , index,Q);
    end
    
else % calculates pulse amplitude for given flip angle for a rectangular pulse
    nu1=beta/(2*pi*1e6*t);
end

end

