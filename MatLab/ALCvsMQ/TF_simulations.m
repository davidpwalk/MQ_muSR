clear system
clear options
clear sequence

%-- Settings --%
save_all_data = false;

% Zeeman (gamma / GHz/T)
ge = -28.02495;
gmu = 0.1355;

% Coupling constants (in GHz) and rotation angles (in degree)
A_iso = 0;
D_parallel = 0.03;
D_perpen = -D_parallel/2;

thetas = deg2rad(linspace(0, 90, 800));
% thetas = deg2rad([1, 5, 20, 45, 70, 85, 89]);
% thetas = deg2rad([45]);
phis = deg2rad([0]); % Phi has no impact on the spectra

% Range of B0
% magnetic_fields = linspace(0, 0.4, 1000);
% magnetic_fields = [1.82, 1.8922];
magnetic_fields = [10];

% the system
system.sqn=[0.5 0.5];       % spin quantum numbers

system.interactions = {};
                     
system.init_state='ex'; % LF mode (muon in the Iz eigenstate)
system.eq = 'init';  % equilibrium state is the same as initial state

% Set options
options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex', 'xe'};
options.labframe = 1;       % lab frame simulation is on
options.awg.s_rate = 5;   % gives sampling rate of simulation in GHz

sequence.tp=8000.0;     % vector with event lengths in ns
sequence.detection=ones(1,length(sequence.tp)); % detection always on

%-- Generation of relevant matrices --%
Sx = sop([1/2 1/2],'xe');
Sy = sop([1/2 1/2],'ye');
Sz = sop([1/2 1/2],'ze');

Ix = sop([1/2 1/2],'ex');
Iy = sop([1/2 1/2],'ey');
Iz = sop([1/2 1/2],'ez');

% Rotation matrices
Rz = @(phi) [cos(phi) -sin(phi) 0;
               sin(phi)  cos(phi) 0;
               0           0          1];

Ry = @(theta) [cos(theta)  0 sin(theta);
             0         1        0;
             -sin(theta) 0 cos(theta)];

T_principal = diag([D_perpen, D_perpen, D_parallel]);

% Number of combinations
Norient = numel(thetas) * numel(phis);

% Preallocate cell array
T_labs = cell(1, Norient);

idx = 1;
for n = 1:length(thetas)
    for m = 1:length(phis)
        theta = thetas(n);
        phi   = phis(m);

        % rotation
        R = Rz(phi) * Ry(theta);

        % rotated tensor
        T_labs{idx} = R * T_principal * R.';

        idx = idx + 1;
    end
end

%% Simulation
[experiment, options] = triple(sequence, options);  % build experiment to get experiment.tp

% Preallocate: one cell per orientation. Each cell holds a matrix: (#det_ops) x Nfields
Nfields = numel(magnetic_fields);
Nt = experiment.tp/experiment.dt + 1;

det_op = 3;  % 1=ez, 2=ex, 3=xe

eigenvalues = cell(1, Norient);
spectra = zeros(Norient, numel(options.det_op), Nfields);
powder_signal = zeros(1, Nt);

weights = sin(thetas);
weights = weights/sum(weights);

if save_all_data
    signals = cell(1, Norient);
    allsignals = cell(Norient, Nfields);     % only if you need full detected_signal per field
end

tic;

% parpool('Threads')  % Start ThreadPool, threads share memory
for k = 1:Norient
    T_lab = T_labs{k};
    temp_eigenvalues = zeros(4, Nfields);

    if save_all_data
        temp_signal = zeros(length(options.det_op), Nt, Nfields);
        temp_allsignals = cell(1, Nfields);
    end

    parfor n = 1:Nfields
        temp_nu_muon = gmu * magnetic_fields(n);
        temp_nu_electron = ge * magnetic_fields(n);

        temp_system = system;  % copy
        temp_system.interactions = { ...
            1,0,'z','e', temp_nu_electron; ...
            2,0,'z','e', temp_nu_muon; ...
        };

        [temp_experiment, temp_options] = triple(sequence, options);
        [temp_system, temp_state, temp_options] = setup(temp_system, temp_options);

        % Isotropic dipolar Hamiltonian
        A_iso_matrix = diag([A_iso, A_iso, A_iso]);
        H_iso = Sx*A_iso_matrix(1, 1)*Ix + Sy*A_iso_matrix(2, 2)*Iy + Sz*A_iso_matrix(3, 3)*Iz;

        % Dipolar Hamiltonian in lab frame
        H_dip = Sx*T_lab(1,1)*Ix + Sx*T_lab(1,2)*Iy + Sx*T_lab(1,3)*Iz + ...
                Sy*T_lab(2,1)*Ix + Sy*T_lab(2,2)*Iy + Sy*T_lab(2,3)*Iz + ...
                Sz*T_lab(3,1)*Ix + Sz*T_lab(3,2)*Iy + Sz*T_lab(3,3)*Iz;

        temp_system.ham = temp_system.ham + 2*pi*(H_dip + H_iso);

        temp_eigenvalues(:, n) = eig(temp_system.ham);

        [temp_state, detected_signal, ~] = homerun(temp_state, temp_system, temp_experiment, temp_options, []);

        if save_all_data
            temp_signal(:, :, n) = detected_signal.sf;
            temp_allsignals{n} = detected_signal;
        end

        powder_signal = powder_signal + weights(k) * detected_signal.sf(det_op, :);
    end
    
    eigenvalues{k} = temp_eigenvalues;
    % spectra(k,:,:) = temp_integrals;

    if save_all_data
        signals{k} = temp_signal;
        allsignals(k,:) = temp_allsignals; 
    end
end

toc;

%% Fourier Transform to get powder spectrum
time_axis = linspace(0, experiment.tp, Nt);
window = cos(pi/2 * (0:length(time_axis)-1) / (length(time_axis)-1));

% Apply apodization to signal
powder_signal_apod = powder_signal .* window;

% figure;
% plot(time_axis, real(powder_signal), 'k', 'DisplayName', 'Original')
% hold on
% plot(time_axis, real(powder_signal_apod), 'r', 'DisplayName', 'Apodized')
% xlabel('Time (\mus)')
% ylabel('Signal (a.u.)')
% legend
% title('Time-Domain Signal with Apodization')

signal_fft = fftshift(fft(powder_signal_apod));
freq_axis = linspace(-options.awg.s_rate/2, options.awg.s_rate/2, Nt); % GHz

fig = figure('NumberTitle','off','Name','TF Powder Spectrum');
subplot(2,1,1)
plot(time_axis, real(powder_signal_apod), 'k')
xlabel('Time / ns')
ylabel('Signal (a.u.)')
title('Weighted Time-Domain Powder Signal (TF)')

subplot(2,1,2)
plot(freq_axis, abs(signal_fft), 'r')
xlabel('Frequency (GHz)')
ylabel('Amplitude (a.u.)')
title('TF Powder Spectrum of Muonium')

%% Plot time evolution of signal
% [experiment, options] = triple(sequence, options);  % build experiment to get experiment.tp

signal = zeros(size(signals{1}));
Nt = size(signals{1}, 2);
time = (0:Nt-1) * experiment.dt;

if save_all_data
    stride = 1;   % downsample
    t_idx = 1:stride:length(time);
    trace = squeeze(signals{1}(1, t_idx, 1));
    time_ds = time(t_idx);
    
    fig = figure('NumberTitle','off','Name','Time-Domain Spectrum Pz');
    hold on
    plot(time_ds, real(trace))
    xlabel('Time / ns')
    ylabel('Polarization')

    % save('Data/ALC_signal_time_evolution_off_resonance.mat', 'trace', 'time_ds')
end
