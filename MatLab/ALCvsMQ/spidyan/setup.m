function [system, state, options]=setup(system, options)
%  [system, state, options] = setup(system, options)
%
%  The function setup processes the input parameters of system, so that 
%  they can be used to call the function homerun. The input must contain a
%  vector with spin quantum numbers .sqn of the present spins. The ordering
%  within this vector will apply to every other option which is spin
%  specific. The Hamiltonian .ham is built from an interaction table (for
%  examples see below). 
%  The initial state can be built from a string, allowing to set the state 
%  of each spin individually, or can be given as matrix The equilibrium 
%  state can be left empty (assumption of high temperature approximation 
%  for  all spins), defined individual through a string or to set to use 
%  the initial staet with the string 'init'.
%  If a spin with S>1/2 is found, the program looks for D, E and angles
%  theta and phi. These can be defined as vectors .D,.E,... (position
%  correspond to positions in .sqn. If vectors have a length of 1, the same
%  D, E, phi and theta are assumed for all high-spins. Alternatively, the
%  ZFS can be constructed with the interaction table.
%  If relaxtion is switched on, the program calculates the relaxtion super
%  operator .gamma from provided .T1 and .T2 times. It is possible to
%  define transition selective relaxation times. For a discussion of that,
%  please refer to the documentation. Non defined transitions are
%  substituted with default values, which can either be set manually as
%  .defaultT1 and .defaultT2 or are assumed by the program as very large.
%  The detection operators are build from strings. Input in matrix form is
%  also possible. Transition selective operators can also be interpreted.
%  See below and documentation for more details.
% 
%  Examples 
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
%  Equlibrium state and initial state are both strings, or alternatively
%  matrices. In the former, the state of each spin can be specified with
%  'z', 'x', ...
% 
%  system.init = 'zz'  % both spins are quantizied along z
%  system.init = 'zx'  % first spin along z, second in coherence
%  system.init = 'ze'  % first spin along z, second not defined, for
%                        example the case in electron nucleus system
%  system.init = 'z'   % same as above in the case of a two (or more) spin
%                        system, or in a  single spin system, only the
%                        first spin is quantizied along z
%  Also possible are:
%  system.init = '-zx'
%  system.init = '-z-x'
%  system.init = 'z-z'
%  system.init = 'ex'
% 
%  The equilibrium state (required for relaxation), can be set in the same
%  manner as the initial state
% 
%  system.eq = '-zz'
% 
%  or can be the same as the initial state:
% 
%  system.eq = 'init'
% 
%  Relaxation times can either be given as a single or in matrix form. 
%  In the former case, the same value is applied all T1 or T2 relaxation
%  pathways:
%  
%  system.T1 = 1000;
%  system.T2 = 500;
%
%  The latter allows to define T1 and T2 values to all transitions 
%  seperately. The used matrix has the same structure as a denisty matrix
%  of the system, with elements corresponding to the same states as in
%  sigma. The off-diagonal elements correlate the pathways between the
%  states in the first dimension and the second dimension. Therefore it is
%  necessary to define the upper triangle only.
% 
%  Assign different relaxation times to all transitions in above example:
% 
%                aa  ab     ba      bb
%  system.T1 =  [0   10^7   10^3    10^4;  aa
%                0   0      10^4    10^3;  ab
%                0   0      0       10^7;  ba
%                0   0      0       0   ;] bb
% 
%                aa  ab     ba      bb
%  system.T2 =  [0   10^3   10^2    10^3;  aa
%                0   0      10^3    10^2;  ab
%                0   0      0       10^5;  ba
%                0   0      0       0   ;] bb
% 
%  From this, the program builds the relaxation superator .gamma in 
%  Liouville space.
% 
%  St. Pribitzer, 2015

% Default relaxation times in ns, used for non-defined transitions. Very 
% large relaxation times, basically minimizing relaxation
default_T1=1000000;
default_T2=1000000;

% checks if for spins are given at all
if ~isfield(system,'sqn') || isempty(system.sqn)    
    error('Setup.m called without spin quantum number.');
end

% determines number of spins
if ~isfield(system,'spins')
    system.spins = length(system.sqn);
end

% checks if for some reason negative sqn are given
for i = 1: length(system.sqn)
    if sign(system.sqn(i))<0
        system.sqn(i)=abs(system.sqn(i));
        warning('Negativ spin quantum number for spin %1.0f. Using %1.1f instead.',i,system.sqn(i));
    end
end

% checks if interaction table is consistent with given spins
for i = 1:size(system.interactions,1)
    if any(cell2mat(system.interactions(i,1:2))>system.spins)
        error('More spins in interaction table than defined in system.sqn.')
    end
end

% calls interactive interface for inputing relaxation times, might be
% dropped in future versions
if isfield(options,'relaxation') && options.relaxation &&  isfield(options,'ui') && ~isempty(options.ui) && options.ui
    [system] = ui_relaxation(system, options);
end

% if no default T1 and T2 are given, the previously defined values are
% assumed
if ~isfield(system,'default_T1') || isempty(system.default_T1)
    system.default_T1=default_T1;
end
if ~isfield(system,'default_T2') || isempty(system.default_T2)
    system.default_T2=default_T2;
end

% if no Hamiltonian is given, this checks if the input is sufficient to
% calculate one
if ~isfield(system,'ham') && ~isfield(system,'interactions')
    error('Neither Hamiltonian nor interaction table available.');
end

% if no Hamiltonian is available, its created here by calling
% build_hamiltonian with the interaction table
if ~isfield(system,'ham') 
    system = build_hamiltonian(system,options);
end

% Creation of initial state from string or matrix, or if none given, high
% the temperature limit is assumed for all spins.
if  ~isfield(system,'init_state')
    system.init_state=[];
end
state = init_state(system.init_state,system,'init');

% creates equilibrium state from string or matrix. If none given, assumes
% that initial state is in equilibrium. Can be set to be the same as
% initial state with 'init'

if ~isfield(system,'eq') || (~isnumeric(system.eq) && strcmp(system.eq,'init')) || isempty(system.eq)
    if ~isfield(system,'eq') || isempty(system.eq)
        fprintf('Setup.m called without equilibrium state, assuming initial state is equilibrium state.\n');
    end
    system.eq = state;
else
    system.eq = init_state(system.eq,system,'eq');
end
    
% If not defined, relaxation is switched off
if ~exist('options','var') || ~isfield(options,'relaxation') || isempty(options.relaxation)
    options.relaxation=0;
    fprintf('Relaxation not defined, switching it off.\n');
end

% creates relaxation superoperator from relaxation times by calling
% changeup
if options.relaxation && ~isfield(system,'gamma')
    if isfield(system,'T1') && isfield(system,'T2') 
            system.gamma=changeup(system);
    else
        error('Either T1 or T2 are missing to compute relaxation superoperator.')  
    end
elseif ~options.relaxation
    system.gamma=[];
end

% compuation of detection matrices from strings. Input can also be a
% matrix. Empty if no detection operators are given
if ~isfield(options,'det_op') || isempty(options.det_op)
    fprintf('No detection operator given.\n');
    options.detect=[];
else
    %creates detection operators from strings
    for k = 1: length(options.det_op)
        options.detect{k}=build_operators(system,options.det_op{k});
    end
end

% looks for highest frequency in the interaction table and stores it. Can
% then be used for automatic down conversion, in case down conversion
% frequency is not explicitly given
system.highest_freq = 0;
for k = 1 : size(system.interactions,1)
    if abs(system.interactions{k,5}) > system.highest_freq
        system.highest_freq = system.interactions{k,5};
    end
end

    
    

% % % % % % % needs editing
% if options.truncation is true, SPIDYAN searches for transitions outside
% the Nyqist range, and sets the corresponding elements in the excitation
% operator to zero.
% if isfield(options,'truncation') && ~isempty(options.truncation)
%     low_cut = find(diff(diag(system.ham))/2/pi>0,1,'first');
%     high_cut = find(diff(diag(system.ham))/2/pi<-options.awg.s_rate/2,1,'last');
%     
%     if (low_cut < high_cut)
%         bla = high_cut;
%         high_cut = low_cut;
%         low_cut = bla;
%     end
%     
%     if ~isempty(low_cut)
%         system.sops.Sx(low_cut:end,low_cut:end) = zeros((2*system.sqn1+1)-low_cut+1);
%     end
%     if ~isempty(high_cut)
%         system.sops.Sx(1:high_cut+1,1:high_cut+1) = zeros(high_cut+1);
%     end
%     
% end

end

function state = init_state(input,system,type)
% Local function init_state which creates the inital state and the
% equilibrium state from strings.

if  isempty(input) 
    fprintf('Setup.m called without initial state, assuming high temperature limit.\n');
    string = '-z';
    state = build_operators(system,string);
    for k = 2: length(system.sqn)
        string(k) = 'e';
        string(k+1) = 'z';
        sigma=build_operators(system,string);
        state = sigma + state;
    end
elseif ~isnumeric(input)
    state = zeros(size(system.ham));
    string = char(length(system.sqn));   
    for k = 1 : length(system.sqn)
        string(k) = 'e';
    end
    fac = 1;
    index = 1;
    for k = 1: length(input)
        if strcmp(input(k),'-')
            fac = -1;
        elseif strcmp(input(k),'+')
            fac = 1;
        elseif strcmp(input(k),'0')
            index = index+1;
        else
            stringl=string;
            stringl(index)=input(k);
            sigma = build_operators(system,stringl);
            state = fac*sigma + state;
            index = 1+index;
            fac = 1;
        end
        if index > length(system.sqn)
            break
        end
    end 
elseif isnumeric(input) && size(input,1) == size(input,2) && size(input,1) ~= 1
    state=input;
else
    switch type
        case 'init'
            error('Initial state not recognized.')
        case 'eq'
            error('Equilibrium state not recognized.')
    end
end

end
