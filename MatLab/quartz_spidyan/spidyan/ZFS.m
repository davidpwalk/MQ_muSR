function [H_tot, H_zfs] = ZFS( S,D,E,theta,phi,H0 )
% ZFS( S,D,E,Dpa,H0 ) constructs the Zerofield Splitting Hamiltonian in the Eigenframe
% of the Zeemann Hamiltonian 
% constructs perturbated Hamiltonian if H0 (zeeman) is given
% 
% INPUT
% S: Spin Quantum Number
% D, E: ZFS Parameters
% theta: arcsin(z/r)
% phi:  atan2(y,x)
% H0: zero order Hamiltonian (e.g. zeeman)
% 
% OUTPUT
% H_zfs: zfs Hamiltonian rotated by theta and phi euler matrix
% H_tot: total perturbated Hamiltonian in Eigenframe of Zeeman
% 
% Andrin Doll and Nino Wili, 2014

H_zfs=zeros(2*S+1);

% build Operators
Sx = spops(S,'x');
Sy = spops(S,'y');
Sz = spops(S,'z');
Sxyz(1) ={Sx} ;
Sxyz(2) ={Sy} ;
Sxyz(3) ={Sz} ;

%construct Dtensor in the Eigenframe of zeeman
Dpv = [-1,-1,2]/3*D + [1,-1,0]*E ;
% rotation matrix 
% Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
% Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1]; 
% Rzprime = Rx'*Rz*Rx; %transformation into ""
R=rotmat(phi,theta,0);
Dtensor = R*diag(Dpv)*R';   % computation of the full D tensor

%calculate (Sx,Sy,Sz)Dtensor(Sx,Sy,Sz)'
for ii = 1:3
   for jj = 1:3 
    H_zfs = H_zfs + Dtensor(ii,jj)*Sxyz{ii}*Sxyz{jj}; 
   end
end

%%%%%%%%%%%%%%
%Perturbation
%%%%%%%%%%%%%%
if(nargin>4)

% get zeroth and first order Energies
for n = 1:(2*S+1)
   
    E0(n) = H0(n,n);
    E1(n) = H_zfs(n,n);
    
end

%get second order energies
E2=zeros(1,(2*S+1));

for n = 1:(2*S+1)
   for m = 1:(2*S+1) 
   if(n~=m)
   
       E2(n) = E2(n) + (H_zfs(m,n)*H_zfs(n,m))/(E0(n)-E0(m));
       
   end
   end
end

%make new Hamiltonian
for n = 1:(2*S+1)
    Energy(n)=E0(n) + E1(n) + E2(n);
%      Energy(n)=E0(n) + E1(n);  %only first order
end

Energy = Energy';

H_tot = diag(Energy);

end
end

