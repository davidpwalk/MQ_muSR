function system = build_hamiltonian(system,options)
% system = build_hamiltonian(system,options)
% Function build_hamiltonian creates the free evolution Hamiltonian for a
% given spin system. The input consits of the two structures system and
% options. The structure system must contain a vectore .sqn with the spin
% quantum numbers and an interaction table .interactions, from which the
% Hamiltonian is built up. The structure options is optional. If .sqn
% contains spins with S > 1/2 ZFS is included through pertubation theory
% and options must provide the fields .labframe and .LO. The former tells 
% the program if it should build the Hamiltonian at the simulation or lab 
% frame frequency. The latter is the frequency used for upconversion from 
% the simulation to the lab frame. 
% 
% Examples 
%  Input:
%  Spin quantum numbers are defined as vector
% 
%  system.sqn = [1st 2nd 3rd ...]
% 
%  for example a two-spin system with two S = 1/2 spins:
% 
%  system.sqn = [1/2 1/2]
% 
%  The interaction table contains all Hamiltonian terms, and has the
%  following structure
% 
%  system.interactions = {index one, index two, operator one, op two, size}
% 
%  The operators are given as strings ('x', 'y', 'z', 'p', 'm', 'e') or as
%  matrix. For example, the Zeeman terms Sz of aboves example would be 
%  written as
% 
%  system.interactions = {1, 0, 'z', 'e', 1; % 1st spin, w = 1 GHz
%                         2, 0, 'z', 'e', 0.015} % 2nd spin, w = 15 MHz
% 
%  Interspin interactions can be added in a similar matter, for example the
%  hyperfine interaction between a nucleus and an elecron spin, A SzIz and 
%  B SzIx are written as:
% 
%  system.interactions = {1, 0, 'z', 'e', 1; 
%                         2, 0, 'z', 'e', 0.015;
%                         1, 2, 'z', 'z', 0.04;  % SzIz term, with A=40 MHz
%                         1, 2, 'z', 'z', 0.02;} % SzIx term, with B=20 MHz
% 
% St. Pribitzer, 2015

interactions = system.interactions;
system.ham = 0;

% If the program does not find options or either options.LO or
% options.labframe are missing, they are replaced by default values (0)
if ~exist('options','var') || isempty(options)
    options.labframe = 0;
    options.LO = 0;
else
    if  ~isfield(options,'labframe') || isempty(options.labframe)
        options.labframe=0;
    end
    if  ~isfield(options,'LO') || isempty(options.LO)
        options.LO=0;
    end
end

% Check if spins with S > 1/2 are present, and if yes, if ZFS parameters
% are availabe and sets wether or not to compute ZFS with pertubation
% theory.
if any(system.sqn>1/2) && isfield(system,'D') && isfield(system,'E') && isfield(system,'theta') && isfield(system,'phi')
    assumption = 'ZFS';
else
    assumption = 'noZFS';
end


% Computation of the Hamiltonian, takes into account pertubation therory
% for ZFS for Zeeman term (case ZFS), or no ZFS, which strictly processes 
% the interaction table only. The ZFS must always be computed at the
% correct frame, which makes it necessary to provide the fequency
% difference between the simulation frame and the lab frame (options.LO).
% For simulations in the labframe one additionaly has to switch on
% options.labframe by setting it to 1. The program the upconverts all
% Zeeman terms to labframe frequencies. If ZFS is requested in the
% simulation frame, the ZFS interactions are downconverted to the
% simulation frame after pertubation theory. No plausibility check (certain
% transitions might end up in the negative frequency range or outside the
% Nyquist band.
switch assumption
    case 'noZFS'
        for k = 1:size(interactions,1)
            op1 = 1;
            op2 = 1;
            
            if options.labframe
                if ~isempty(system.interactions{k,3}) && strcmp(system.interactions{k,3},'z') && (isempty(system.interactions{k,4}) || strcmp(system.interactions{k,4},'e')) && system.interactions{k,5}>0.1
                    interactions{k,5}=interactions{k,5}+options.LO;
                    fprintf('Adapted interaction %1.0f to labframe.\n',k)
                elseif ~isempty(system.interactions{k,4}) && strcmp(system.interactions{k,4},'z') && (isempty(system.interactions{k,3}) || strcmp(system.interactions{k,3},'e')) && system.interactions{k,5}>0.1
                    interactions{k,5}=interactions{k,5}+options.LO;
                    fprintf('Adapted interaction %1.0f to labframe.\n',k)
                end
            end
                        
            for l = 1:length(system.sqn)
                if isempty(interactions(k,1)) || interactions{k,1}==0 || interactions{k,1}~=l
                    op1 = kron(op1,spops(system.sqn(l),'e'));
                elseif interactions{k,1} == l;
                    op1 = kron(op1,spops(system.sqn(l),interactions{k,3}));
                end
            end
            
            for l = 1:length(system.sqn)
                if isempty(interactions(k,2)) || interactions{k,2}==0 || interactions{k,2}~=l
                    op2 = kron(op2,spops(system.sqn(l),'e'));
                elseif interactions{k,2} == l;
                    op2 = kron(op2,spops(system.sqn(l),interactions{k,4}));
                end
            end
            
            system.H{k} = interactions{k,5}*op1*op2;
            system.ham = 2*pi*system.H{k}+system.ham;
        end
        
    case 'ZFS'
        
        if numel(system.D)<sum(system.sqn>1/2)
            warning('Number of spins with S>1/2 exceeds available ZFS parameters. Assuming the same ZFS for all spins.')
            system.D(1:numel(system.sqn)) = system.D;
            system.E(1:numel(system.sqn)) = system.E;
            system.theta(1:numel(system.sqn)) = system.theta;
            system.phi(1:numel(system.sqn)) = system.phi;
        end
        
        for k = 1:size(interactions,1)
            op1 = 1;
            op2 = 1;
            zfs = 0;
            
            if options.labframe
                if ~isempty(system.interactions{k,3}) && strcmp(system.interactions{k,3},'z') && (isempty(system.interactions{k,4}) || strcmp(system.interactions{k,4},'e')) && system.interactions{k,5}>0.1
                    interactions{k,5}=interactions{k,5}+options.LO;
                    fprintf('Adapted interaction %1.0f to labframe.\n',k)
                    spinpos = 1;
                elseif ~isempty(system.interactions{k,4}) && strcmp(system.interactions{k,4},'z') && (isempty(system.interactions{k,3}) || strcmp(system.interactions{k,3},'e')) && system.interactions{k,5}>0.1
                    interactions{k,5}=interactions{k,5}+options.LO;
                    fprintf('Adapted interaction %1.0f to labframe.\n',k)
                    spinpos = 2;
                end
                if system.sqn(interactions{k,spinpos})>1/2
                    zfs = 1;
                    dwnconv = 0;
                end
            elseif options.LO
                if ~isempty(system.interactions{k,3}) && strcmp(system.interactions{k,3},'z') && (isempty(system.interactions{k,4}) || strcmp(system.interactions{k,4},'e')) && system.interactions{k,5}>0.1
                    interactions{k,5}=interactions{k,5}+options.LO;
                    spinpos = 1;
                elseif ~isempty(system.interactions{k,4}) && strcmp(system.interactions{k,4},'z') && (isempty(system.interactions{k,3}) || strcmp(system.interactions{k,3},'e')) && system.interactions{k,5}>0.1
                    interactions{k,5}=interactions{k,5}+options.LO;
                    spinpos = 2;
                end
                if system.sqn(interactions{k,spinpos})>1/2
                    zfs = 1;
                    dwnconv = 1;
                    fprintf('Calculating ZFS (if available) of interaction %1.0f in labframe.\n',k)
                end
            end
                
            
            if zfs
                
                for l = 1:length(system.sqn)
                    if interactions{k,spinpos}~=l
                        op1 = kron(op1,spops(system.sqn(l),'e'));
                    elseif interactions{k,spinpos} == l;
                        hamZ = interactions{k,5}*spops(system.sqn(l),'z');
                        
                        hamtot=ZFS(system.sqn(l),system.D(l),system.E(l),system.theta(l),system.phi(l),hamZ);
                                    
                        if dwnconv
                            hamtot = hamtot - options.LO*spops(system.sqn(l),'z');
                        else
                            
                        end
                        
                        op1 = kron(op1,hamtot);
                    end
                end
                
                system.H{k} = op1;
                
            else
                
                for l = 1:length(system.sqn)
                    if isempty(interactions(k,1)) || interactions{k,1}==0 || interactions{k,1}~=l
                        op1 = kron(op1,spops(system.sqn(l),'e'));
                    elseif interactions{k,1} == l;
                        op1 = kron(op1,spops(system.sqn(l),interactions{k,3}));
                    end
                end
                
                for l = 1:length(system.sqn)
                    if isempty(interactions(k,2)) || interactions{k,2}==0 || interactions{k,2}~=l
                        op2 = kron(op2,spops(system.sqn(l),'e'));
                    elseif interactions{k,2} == l;
                        op2 = kron(op2,spops(system.sqn(l),interactions{k,4}));
                    end
                end
                
                system.H{k} = interactions{k,5}*op1*op2;
                
            end
            
            system.ham = 2*pi*system.H{k}+system.ham;
        end
        
end