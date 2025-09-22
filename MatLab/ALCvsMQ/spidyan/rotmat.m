function Rp=rotmat(alpha,beta,gamma)
%
% function Rp=rotmat(alpha,beta,gamma)
%
% Compute rotation matrix
%
% alpha, beta, gamma    Euler angles in radians
%                       if there is only a single argument alpha, it is 
%                       interpreted as a vector of three Euler angles
% Rp                    rotation matrix
% 
% 2015, Gunnar Jeschke

if nargin == 1,
    gamma = alpha(3);
    beta = alpha(2);
    alpha = alpha(1);
end;
ca=cos(alpha);
sa=sin(alpha);
cb=cos(beta);
sb=sin(beta);
cg=cos(gamma);
sg=sin(gamma);

Rp=[ca*cb*cg-sa*sg,sa*cb*cg+ca*sg,-sb*cg;-ca*cb*sg-sa*cg,-sa*cb*sg+ca*cg,sb*sg;ca*sb,sa*sb,cb];
end
