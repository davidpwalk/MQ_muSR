clear system
clear options
clear sequence

% Zeeman (gamma / MHz/T)
ge = 28024.95;
gmu = -135.5;

% Coupling constants (in MHz)
A_iso = 514.8;
D_parallel = 2;
D_perpen = -D_parallel/2;

% Range of B0
magnetic_fields = linspace(1.8, 2, 100);

Iz = sop([1/2 1/2],'ez');
Sz = sop([1/2 1/2],'ze');

Ix = sop([1/2 1/2],'ex');
Sx = sop([1/2 1/2],'xe');

%%
% the system
system.sqn=[0.5 0.5];       % spin quantum numbers
system.interactions={1,0, 'z','e', nu_electron;       % Zeeman term of electron in simulation frame
                     2,0, 'z','e', nu_muon;           % Zeeman of muon
                     1,2, 'x','x', A_iso;             % isotropic HF
                     1,2, 'y','y', A_iso;             % isotropic HF
                     1,2, 'z','z', A_iso};            % isotropic HF
% TODO: add axially symmetric Dipolar interaction.
                     
system.init_state='ez'; % LF mode (muon in the Iz eigenstate)
system.eq = 'init';  % equilibrium state is the same as initial state

% the options
options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex'}

options.labframe = 1;     % lab frame simulation is on

[system, state, options]=setup(system,options); % build system 

parfor n=1:length(magnetic_fields)
    nu_muon = gmu*magnetic_fields[n];
    nu_electron = ge*magnetic_fields[n];

    
end