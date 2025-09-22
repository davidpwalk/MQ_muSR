function sim = simexp( v, x )
%SIMEXP simulation of an exponential
%    v  input parameter vector
%       v(1) decay constant (rate of decay)
%       v(2) offset, if length(v)==1, an offset of zero is assumed
%    x  x axis, v(1) relates to these values
%    y  data, only real part is fitted
%
% rmsd  root mean square deviation
% sim   simulated data
%
% G. Jeschke, 2012

sim=exp(-v(1)*x)+v(2);

end

