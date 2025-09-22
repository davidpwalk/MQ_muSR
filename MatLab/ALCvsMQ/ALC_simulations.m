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
nu_muon = gmu*magnetic_fields(1);
nu_electron = ge*magnetic_fields(1);

% the system
system.sqn=[0.5 0.5];       % spin quantum numbers
system.interactions={1,0, 'z','e', nu_electron;       % Zeeman term of electron in simulation frame
                     2,0, 'z','e', nu_muon;           % Zeeman of muon
                     1,2, 'x','x', A_iso;             % isotropic HF
                     1,2, 'y','y', A_iso;             % isotropic HF
                     1,2, 'z','z', A_iso;             % isotropic HF
                     1,2, 'x','x', D_perpen;          % axially symmetric HF
                     1,2, 'y','y', D_perpen;          % axially symmetric HF
                     1,2, 'z','z', D_parallel};       % axially symmetric HF
                     
system.init_state='ez'; % LF mode (muon in the Iz eigenstate)
system.eq = 'init';  % equilibrium state is the same as initial state

% the options
options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex'};

options.labframe = 1;     % lab frame simulation is on

experiment.pulse = {};   % no pulses
experiment.tp = linspace(0,10,1000);
experiment.dt = experiment.tp(2) - experiment.tp(1)

[system, state, options]=setup(system,options); % build system 

[state, detected_signal, experiment] = homerun(state, system, experiment, options, []);


%%
% the system
system.sqn=[0.5 0.5];       % spin quantum numbers

system.interactions = {};
                     
system.init_state='ez'; % LF mode (muon in the Iz eigenstate)
system.eq = 'init';  % equilibrium state is the same as initial state

% the options
options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex'};

options.labframe = 1;     % lab frame simulation is on

experiment.pulse = {};   % no pulses
experiment.tp = linspace(0,10,1000);  

% Create cells for storing results
signals=cell(1,length(magnetic_fields));
allsignals=cell(1,length(magnetic_fields));

tic; % Start stopwatch

parfor n=1:length(magnetic_fields)
    temp_nu_muon = gmu*magnetic_fields(n);
    temp_nu_electron = ge*magnetic_fields(n);

    temp_system = system;
    temp_system.interactions={1,0, 'z','e', temp_nu_electron;       % Zeeman term of electron in simulation frame
                     2,0, 'z','e', temp_nu_muon;           % Zeeman of muon
                     1,2, 'x','x', A_iso;             % isotropic HF
                     1,2, 'y','y', A_iso;             % isotropic HF
                     1,2, 'z','z', A_iso;             % isotropic HF
                     1,2, 'x','x', D_perpen;          % axially symmetric HF
                     1,2, 'y','y', D_perpen;          % axially symmetric HF
                     1,2, 'z','z', D_parallel};       % axially symmetric HF

    [temp_system, temp_state, temp_options]=setup(temp_system,options); % build system 
    
    [temp_state, detected_signal, temp_experiment] = homerun(temp_state, temp_system, experiment, temp_options, []);  % propagate for each magnetic field
    
    signals{n} = detected_signal.sf;    % keep signal from simulation frame

    allsignals{n} = detected_signal;
end

toc;

signal = zeros(size(signals{1}));
time = experiment.tp;

% Integrate
integrals = zeros(length(options.det_op), length(magnetic_fields));

% Filter
filtered_signals = cell(1, length(magnetic_fields));

filter_cutoff = 1.0;    % GHz

cutoff = filter_cutoff/((1/(time(2)-time(1)))/2);

% Filter is a lowpass butterworth filter of 3rd order
N=3;
[B,A] = butter(N,cutoff,'low');

for n=1:length(magnetic_fields)
    signal = signal + signals{n};

    integrals(:, n) = sum(real(signals{n}), 2)/length(signal);

    filtered_signals{n} = signals{n};

    for m=1:size(signals{n}, 1)
        N_pad = 100;
        temp_signal = filtered_signals{n}(m, :);
        len = length(temp_signal);
        temp_signal = [temp_signal(N_pad+1:-1:1+1), temp_signal, temp_signal(len-N_pad+1-1:len-1)];
        temp_filtered_signal = filtfilt(B, A, temp_signal);
        filtered_signals{n}(m, :) = temp_filtered_signal(N_pad+1:len+N_pad);
    end
end
