function [rmsd, sim] = rmsexp( v, x, y )
% RMSEXP root mean square deviation (and optionally) fit of an exponential
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

sim = simexp( v, x );

c = sum(sim.*y)/sum(sim.*sim);

sim=c*sim;
diff=y-sim;
rmsd=sqrt(sum(diff.^2)/length(diff));

end

