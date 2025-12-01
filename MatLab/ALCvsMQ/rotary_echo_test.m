clear system
clear options
clear sequence

% Example from spidyan 2 manual
c=1.4; % center frequency of Gaussian distribution, in GHz
g=0.01; % width of the Gaussian distribution in GHz
nu0_s=1.3; % starting value for sampling of Gaussian distribution
nu0_e=1.5; % final value for sampling of Gaussian distribution
st=0.001; % stepsize for sampling of Gaussian distribution
nu0_vec=nu0_s:st:nu0_e; % vector with resonance frequencies
p=exp(-((c-nu0_vec)/g).^2); % probabilities for each spin
p=p/trapz(p); % normalization

% the system
system.sqn=[0.5 0.5]; % spin quantum numbers
system.interactions={1,0, 'z','e',1.4; % Zeeman term of electron, value is replaced in loop
                     2,0,'z','e',-0.014; % Zeeman of nucleus
                     1,2, 'z','z',-0.045; % secular part of hyperfine interaction , A*SzIz
                     1,2, 'z','x',0.008}; % pseudo secular part of hyperfine interaction , B*SzIz

system.init_state='-ze'; % only the electron spin is quantizised along z, nucleus is not defined, same statement would be '-z'
system.eq = 'init'; % same as initial state

% the options
options.relaxation=0; % no relaxation
options.down_conversion=0; % downconversion of signal off
options.det_op={'pe'}; % detection operator

% the sequence
sequence.tp=[50,300,25,600]; % vector with event lengths in ns
sequence.beta(1)=pi/2; % flip angle for first event
sequence.beta(3)=pi; % flip angle for third event
sequence.frq={[1.3 1.5]}; % initial and final frequency
sequence.detection=ones(1,length(sequence.tp)); % detects all events
sequence.t_rise=10; % rise time of chirp pulses
sequence.excite={[1 0]}; % only excite the first spin
[experiment,options] = triple(sequence, options);

signalc=cell(1,length(nu0_vec)); % creates an empty cell to store results
parfor k=1:length(nu0_vec)
    systeml=system; % copy system into loop
    systeml.interactions{1,end}= nu0_vec(k); % add loop resonance frequency
    [systeml, state, optionsl]=setup(systeml,options); % build system
    [state, signal, experimentl]=homerun(state,systeml,experiment, optionsl,[]); % propagate
    signalc{k}=signal.sf; % keep signals from simulation frame
end

signal=zeros(size(signalc{1})); % creates empty signal

for k=1:length(nu0_vec)
    signal=signal+signalc{k}*p(k); % combines all signals
end
t=linspace(0,sum(experiment.tp),length(signal)); % creates time axis

figure(1); clf; plot(t,real(signal))
title('Signal in simulation frame')
xlabel('time [ns]')
ylabel('<S f1pg>')

signaldc=strike(signal, t, c, options); % down conversion of merged signals
signaldc=signaldc/max(max(signaldc));

figure(2); clf;
title('Signal after down conversion')
hold on
plot(t,real(signaldc))
xlabel('time [ns]')
ylabel('<S f1pg>')
plot([experiment.tp(1) experiment.tp(1)],[-1 1],'r--')
plot([experiment.tp(2)+experiment.tp(1) ...
experiment.tp(2)+experiment.tp(1)],[-1 1],'r--')
plot([experiment.tp(3)+experiment.tp(2)+experiment.tp(1)+experiment.tp(3)+experiment.tp(1)+experiment.tp(2)],[-1 1],'r--')







%% Apply to my problem
clear system
clear options
clear sequence

% Zeeman (gamma / GHz/T)
ge = -28.02495;
gmu = 0.1355;

magnetic_field = 0.0957;

nu_muon = gmu*magnetic_field;
nu_electron = gme*magnetic_field;

% Introduce anisotropy in B1 (change sequence.nu1 in loop)


nu_uw = abs(nu_electron);

figure('NumberTitle','off','Name','nu_electron Anisotropy'); clf;
plot(nu_electron_vec, nu_electron_probs)
xlabel('\nu_e / GHz')
ylabel('Probablility')

% Coupling constants (in GHz) and rotation angles (in degree)
A_iso = 0.003;
D_parallel = 0.0155;
D_perpen = -D_parallel/2;

theta = deg2rad(45);
phis = deg2rad(0);

magnetic_field = 0.0957;

% the system
system.sqn=[0.5 0.5];       % spin quantum numbers

system.interactions = {};
                     
system.init_state='ez'; % LF mode (muon in the Iz eigenstate)
% Rotary only really work in x and y, but one stops the decay of the signal
system.eq = 'init';  % equilibrium state is the same as initial state

% Set options
options.relaxation=0;       % should be off here because decay of signal should be due to inhomogenties
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez'};
options.labframe = 1;       % lab frame simulation is on
options.awg.s_rate = 6;   % gives sampling rate of simulation in GHz

% Set sequence
tau = 8000;  % ns, length of each pulse block

nu_uw = 0.0002;

sequence.tp   = [tau, tau];             % two equal pulse blocks
sequence.nu1  = [1.0, 1.0];             % same amplitudes
sequence.phase = [0, pi];               % phase shift for second part of pulse
sequence.frq  = [nu_uw, nu_uw];         % same frequencies
sequence.t_rise = [0, 0];               % no chirp
sequence.detection = ones(1,length(sequence.tp));

plot_pulses(triple(sequence, options), options)

%-- Generation of relevant matrices --%
Sx = sop([1/2 1/2],'xe');
Sy = sop([1/2 1/2],'ye');
Sz = sop([1/2 1/2],'ze');

Ix = sop([1/2 1/2],'ex');
Iy = sop([1/2 1/2],'ey');
Iz = sop([1/2 1/2],'ez');

Ry = @(theta) [cos(theta)  0 sin(theta);
             0         1        0;
             -sin(theta) 0 cos(theta)];

T_principal = diag([D_perpen, D_perpen, D_parallel]);

T_lab = Ry(theta) * T_principal * Ry(theta).';

% Setup Interactions
system.interactions = { ...
        1,0,'z','e', nu_electron_center; ... % this will be changed in simulation loop anyway
        2,0,'z','e', nu_muon; ...
    };

H_iso = Sx*A_iso*Ix + Sy*A_iso*Iy + Sz*A_iso*Iz;

H_dip = Sx*T_lab(1,1)*Ix + Sx*T_lab(1,2)*Iy + Sx*T_lab(1,3)*Iz + ...
        Sy*T_lab(2,1)*Ix + Sy*T_lab(2,2)*Iy + Sy*T_lab(2,3)*Iz + ...
        Sz*T_lab(3,1)*Ix + Sz*T_lab(3,2)*Iy + Sz*T_lab(3,3)*Iz;

[experiment, options] = triple(sequence, options);
[system, ~, ~] = setup(system, options); % needs to be done so system.ham is defined

system.ham = system.ham + 2*pi*(H_dip + H_iso);



%% Simulation 
signals=cell(1,length(nu_electron_vec)); % creates an empty cell to store results

tic
parfor ii = 1:length(nu_electron_vec)
    temp_system = system;
    temp_system.interactions{1, end} = nu1_vec(ii);

    [temp_system, temp_state, temp_options] = setup(temp_system, options);

    [temp_state, detected_signal, ~] = homerun(temp_state, temp_system, experiment, temp_options, []);

    signals{ii} = detected_signal.sf; % keep signals from simulation frame
end
toc

%% Merge results and plot them

signal=zeros(size(signals{1})); % creates empty signal

for ii=1:length(nu_electron_vec)
    signal=signal+signals{ii}*p(ii); % combines all signals
end
t=linspace(0,sum(experiment.tp),length(signal)); % creates time axis

figure(1); clf; plot(t,real(signal))
title('Signal in simulation frame')
xlabel('time [ns]')
ylabel('<O>')

% signaldc=strike(signal, t, c, options); % down conversion of merged signals
% signaldc=signaldc/max(max(signaldc));
% 
% figure(2); clf;
% title('Signal after down conversion')
% hold on
% plot(t,real(signaldc))
% xlabel('time [ns]')
% ylabel('<S f1pg>')
% plot([experiment.tp(1) experiment.tp(1)],[-1 1],'r--')
% plot([experiment.tp(2)+experiment.tp(1) ...
% experiment.tp(2)+experiment.tp(1)],[-1 1],'r--')
% plot([experiment.tp(3)+experiment.tp(2)+experiment.tp(1)+experiment.tp(3)+experiment.tp(1)+experiment.tp(2)],[-1 1],'r--')


