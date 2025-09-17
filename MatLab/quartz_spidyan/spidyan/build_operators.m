function [sigma] = build_operators(system,string)
% [sigma] = build_operators(system,string)
% 
% The function build_operators builds operators for arbitrary spin system
% from strings. If the input already is a matrix, the input is returned as
% output. The function allows the user to build a vast variety of multi
% spin operators as well as transition selective operators. 
% The string can contain any combination of x, y, z, p, m, e, + and -
% (and numbers for transition selection). The maximum number of letters
% interpreted by the program is the length of system.sqn, which contains
% the spin quantum numbers of all spins in the system and is required. The
% position of the letters in the string corresponds to the position in
% system.sqn (+ and - do not count). The letters represent Pauli matrices,
% raising (p) and lowering (m) lowering operators and unity (e) matrices.
% The signs + and - are interpreted as such, meaning, that the matrix, of
% the next following letter will be multiplied with 1 or -1. It is not
% necessary to provide a letter for all spins in the system, omitted spins
% are interpreted as e. Transition selective operators are controlled with
% numbers, where the number corresponds to the position in the density
% matrix. If a letter is followed by one number #, this is interpreted as a
% population the letter ignored, and the diagonal element at position # is 
% returned. On the other hand, if the program finds two numbers ##, the 
% preceding letter is taken into account, and the returned density matrix
% has the form of the letter at the position # #. See examples below and 
% the online documentation for more details.
% 
% Examples:
% 
% Assuming a single spin with S = 1/2 the Pauli matrices and raising and
% lowering operators are obtained in the following way:
% 
% system.sqn = 1/2
% string = 'z'      % Sz
% string = 'x'      % Sx
% string = 'y'      % Sy
% string = 'e'      % E
% string = 'p'      % S+
% string = 'm'      % S-
% 
% sigma = build_operators(system,string)
% 
% Lets have a look at a high-spin system with S = 1
% 
% system.sqn = 1
% string = 'z'      % Sz
% 
% Sz now has the form 
%          1   0   0
% sigma =  0   0   0 
%          0   0  -1
% 
% Electron spins are usually detected via:
% 
% string = '-z'
% 
%         -1   0   0
% sigma =  0   0   0 
%          0   0   1
% 
% Lets say we want to detect the population of the third level. In that
% case we have to adapt the string to:
% 
% string = 'z3'
% 
% which returns 
%          0   0   0
% sigma =  0   0   0 
%          0   0   1
% 
% Polarization between the second and third level can be detected with 
% 
% string = 'z23'
% 
% which returns 
%          0   0   0
% sigma =  0  1/2  0 
%          0   0 -1/2
% 
% The following expressions are equal:
% 
% string = 'z23'
% string = '-z32'
% 
% Coherence between two state can be detected with
% 
% string = 'x23'
% 
%          0   0   0
% sigma =  0   0  1/2
%          0   1/2 0
% 
% The same rules of apply to raising and lowering operators:
% 
% string = 'p23'
% 
%          0   0   0
% sigma =  0   0   1
%          0   0   0
% 
% Now lets have a look at a multi-spin system:
% 
% system.sqn = [1/2 1/2 1/2]
% 
% The same rules as above are valid here as well.
% 
% For example S1z2z3z is obtained with
% 
% string = 'zzz'
% 
% S1z:
% 
% string = 'zee'
% 
% S2z:
% 
% string = 'eze'
% 
% S3z:
% 
% string = 'eez'
% 
% And transition selective operators, several possibilites for the same
% operator:
%  
% string = 'xee23'
% 
% string = 'exe23'
% 
% string = 'x23'
% 
% All of the above just detect <Sx> between the second and third element of
% the denisity matrix.
% 
% For simplicity it is also possible to omit letters, for example 
% 
% string = 'zee'
% 
% is the same as
% 
% string = 'z'
% 
% or 
% 
% string = 'eze'
% 
% equals
% 
% string = 'ez'
% 
% If more letters than spin quantum numbers are given, surplus letters are
% ignored:
% 
% string = 'zzzz'
% 
% equals
% 
% string = 'zzz'
% 
% A string can contain multiple signs:
% 
% string = '-z+z-z'
% 
% St. Pribitzer, 2015

if ~isnumeric(string)
sign = 1;
pos = 1;
op = 1;
number = false;
type = 'z';

for k = 1 : length(string)
    if strcmp(string(k),'-')
        sign = -1;
    elseif strcmp(string(k),'+')
        sign = 1;
    elseif pos <= length(system.sqn) && (isnan(str2double(string(k))))
        mat = sign * spops(system.sqn(pos),string(k));
        op = kron(op,mat);
        pos = pos + 1;
        sign = 1;
    elseif str2double(string(k))==1i
        k = k + 1;
        mat = sign * spops(system.sqn(pos),'e');
        pos = pos + 1;
        op = kron(op,mat);
        rest = (string(k:end));
        number = true;
        break
    elseif isnumeric(str2double(string(k)))
        rest = (string(k:end));
        number = true;
        break
    end
end

if pos <= length(system.sqn)
    for p = pos : length(system.sqn)
        mat = spops(system.sqn(p),'e');
        op = kron(op,mat);
        pos = pos + 1;
    end
end

if number
    for j = 1 : k - 1
        if strcmp(string(j),'-') 
            sign = -1;
        elseif ~strcmp(string(j),'z') && ~strcmp(string(j),'+') 
            if strcmp(string(j),'p')
                type = 'p';
            elseif strcmp(string(j),'m')
                type = 'm';
            elseif strcmp(string(j),'x')
                type = 'x';
            elseif strcmp(string(j),'y')
                type = 'y';
            elseif strcmp(string(j),'e')
                type = 'e';
            elseif strcmp(string(j),'i')
                type = 'population';
            end
            
            break
        end
    end
    
    opnew = zeros(size(op));
    pop = 1;
    
    for r = 1 : length(rest)
        
        if strcmp(rest(r),',')
            o = str2double(rest(1:r-1));
            u = str2double(rest(r+1:end));
            pop = 0;
            break
        end
    end
    
    if pop
        o = str2double(rest);
    end
    
    
    switch type
        case 'z'
            opnew(o,o)=1/2;
            opnew(u,u)=-1/2;
        case 'p'
            opnew(o,u)=1;
        case 'm'
            opnew(u,o)=1;
        case 'x'
            opnew(o,u)=1/2;
            opnew(u,o)=1/2;
        case 'e'
            opnew(o,o)=1;
            opnew(u,u)=1;
        case 'y'
            opnew(o,u)=-1/2*1i;
            opnew(u,o)=1/2*1i;
        case 'population'
            opnew(o,o)=1;
    end
            
    op = opnew* sign;
end

sigma = op;

elseif isnumeric(string) && size(string,1)==size(string,2)
    
    sigma = string;
    
else
    error('Detection operator not recognized.')
end

end

