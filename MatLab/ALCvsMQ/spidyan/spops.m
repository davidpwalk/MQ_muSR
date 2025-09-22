function [ma] = spops(spinsystem, string)
% [ma] = spops(spinsystem, string)
% 
% This function creates the spin operator matrix from a string for a given spin 
% system. 
% 
% Input: 
% spinsystem    Spin system, a vector with spin quantum numbers.
% string        A string containing the operatortype for each component of 
%               the spin system. Can be 'x', 'y', 'z', 'e' for identity,    
%               'p' for raising and 'm' for lowering operator, or any 
%               combination for spin systems with more components.
% 
% Output:
% ma            spin operator in matrix form.
% 
% Examples: 
% 
% sops([1/2],'z');       % returns Sz operator for a S=1/2 system
% 
% sops([1/2 1/2],'zx')   % returns SzIx operator for a S=I=1/2 system
% 
% St. Pribitzer, 2014


if strcmp(string(1),'-') % checks sign of the string
    op=string(2:end);
    minus=1;
else
    op=string;
    minus=0;
end

l=length(op);

% checks if the number of spins fits the operator
if length(spinsystem)~=l 
    error('Number of components does not match.');
else
    nsv=spinsystem;
end

spmas=cell(1,l);

for k=1:l 
    % creates the spin matrices for requested type, if string is not
    % recognized, an error is returned
    sn=2*nsv(k)+1;
    if strcmp(op(k),'z')
        opma=zeros(sn,sn);
        e=nsv(k);
        for kk=1:sn
            opma(kk,kk)=e;
            e=e-1;
        end
    elseif strcmp(op(k),'e')
        opma=zeros(sn,sn);
        for kk=1:sn
            opma(kk,kk)=1;
        end
    elseif strcmp(op(k),'x')
        opma=zeros(sn,sn);
        for kk=1:(sn-1)
            e=sqrt(kk*(sn-kk));
            opma(kk,kk+1)=e;
            opma(kk+1,kk)=e;
        end
        opma=opma/2;
     elseif strcmp(op(k),'y')
        opma=zeros(sn,sn);
        for kk=1:(sn-1)
            e=sqrt(kk*(sn-kk));
            opma(kk,kk+1)=-1i*e;
            opma(kk+1,kk)=1i*e;
        end
        opma=opma/2;
    elseif strcmp(op(k),'p')
        opma=zeros(sn,sn);
        for kk=1:(sn-1)
            e=sqrt(kk*(sn-kk));
            opma(kk,kk+1)=e;
        end        
    elseif strcmp(op(k),'m')
        opma=zeros(sn,sn);
        for kk=1:(sn-1)
            e=sqrt(kk*(sn-kk));
            opma(kk+1,kk)=e;
        end
    else 
       error('String not recognized, no spin matrix created.');
    end
    spmas{k}=opma;
end

% combines all the single spin operators to the multispin operator by using
% kronecker deltas
if l==1
    ma=opma;
else
    matrix=spmas{1};
    for k=2:l
        matrix=kron(matrix,spmas{k});
    end
    ma=matrix;
end
        
% setting of sign
if minus
    ma=-ma;
end

end

