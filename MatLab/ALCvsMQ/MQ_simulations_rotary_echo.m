clear system
clear options
clear sequence

%-- Settings --%
rot_echo = true;

save_all_data = true;
save_traces_to_disk = false;  % only applies if save_all_data=true
save_trace_dir = 'Data/traces/test';

if save_traces_to_disk && ~exist("save_trace_dir", 'dir')
    mkdir(save_trace_dir)
end

% Zeeman (gamma / GHz/T)
ge = -28.02495;
gmu = 0.1355;

% Coupling constants (in GHz) and rotation angles (in degree)
A_iso = 0.0;
D_parallel = 0.0155;
% D_parallel = 0;
D_perpen = -D_parallel/2;

thetas = deg2rad(linspace(0, 90, 100));
% thetas = deg2rad([1, 5, 20, 45, 70, 85, 89]);
% thetas = deg2rad([45, 60]);
phis = deg2rad([0]); % Phi has no impact on the spectra

magnetic_field = 0.08178;

nu_muon = gmu * magnetic_field;
nu_electron = ge * magnetic_field;

system.interactions = { ...
    1,0,'z','e', nu_electron; ...
    2,0,'z','e', nu_muon; ...
};

% the system
system.sqn=[0.5 0.5];       % spin quantum numbers
                     
system.init_state='ez'; % LF mode (muon spin polarized)
system.eq = 'init';  % equilibrium state is the same as initial state

% Relaxtion matrix (e- transitions 2 us, rest 1 s)
system.T1 = [0, 1e9,  1e9,  1e9;
             0,   0,  1e9,  1e9;            
             0,   0,    0,  1e9;
             0,   0,    0,    0];

system.T2 = [0, 1e9, 2000, 2000;
             0,   0, 2000, 2000;            
             0,   0,    0,  1e9;
             0,   0,    0,    0];

% Set options
options.relaxation=1;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex', 'ze'};
options.labframe = 1;       % lab frame simulation is on
options.awg.s_rate = 6;   % gives sampling rate of simulation in GHz

%

nu_muon = gmu*B0;
nu_electron = ge*B0;

nu_uw = abs(nu_electron);
% nu_uw = 3000

% Set sequence
tau_array = [50, 100, 200, 400, 800];  % ns

sequence.nu1  = [1.0, 1.0];
% sequence.nu1  = [0, 0];
sequence.frq  = [nu_uw, nu_uw];
sequence.t_rise = [0, 0];               % no chirp

if rot_echo
    sequence.phase = [0, pi];               % phase shift for second part of pulse
else
    sequence.phase = [0, 0];               % no phase shift
end

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
% Preallocate: one cell per orientation. Each cell holds a matrix: (#det_ops) x Nfields

eigenvalues = cell(1, Norient);
spectra = zeros(Norient, numel(options.det_op));

tic;

spectra_taus = zeros(Norient, length(tau_array));
LF_asymmetries = zeros(Norient, length(tau_array));

powder_signals = cell(1, length(tau_array));   % or Nt x Ntau matrix if preferred

% parpool('Threads');  % Start ThreadPool, threads share memory
for tau_idx = 1:length(tau_array)
    tau = tau_array(tau_idx);
    
    sequence.tp   = [tau, tau];             % two equal pulse blocks
    sequence.detection = ones(1,length(sequence.tp));
    
    [experiment, options] = triple(sequence, options);  % build experiment to get experiment.tp
    
    Nt = sum(experiment.tp)/sum(experiment.dt) + 1;

    if save_all_data
        signals = cell(1, Norient);   % will hold ndet x Nt arrays
    end

    parfor n = 1:Norient
        T_lab = T_labs{n};
        temp_eigenvalues = zeros(4, Nfields);
        temp_integrals = zeros(numel(options.det_op), Nfields);
    
        if save_all_data
            temp_signal = zeros(length(options.det_op), Nt, Norient);
            % temp_allsignals = cell(1);
        end
    
        [temp_experiment, temp_options] = triple(sequence, options);
        [temp_system, temp_state, temp_options] = setup(system, temp_options);
    
        % Isotropic dipolar Hamiltonian
        A_iso_matrix = diag([A_iso, A_iso, A_iso]);
        H_iso = Sx*A_iso_matrix(1, 1)*Ix + Sy*A_iso_matrix(2, 2)*Iy + Sz*A_iso_matrix(3, 3)*Iz;
    
        % Dipolar Hamiltonian in lab frame
        H_dip = Sx*T_lab(1,1)*Ix + Sx*T_lab(1,2)*Iy + Sx*T_lab(1,3)*Iz + ...
                Sy*T_lab(2,1)*Ix + Sy*T_lab(2,2)*Iy + Sy*T_lab(2,3)*Iz + ...
                Sz*T_lab(3,1)*Ix + Sz*T_lab(3,2)*Iy + Sz*T_lab(3,3)*Iz;
    
        temp_system.ham = temp_system.ham + 2*pi*(H_dip + H_iso);  % 2pi needed cos they are in angular frequency units in spidyan
        
        [temp_state, detected_signal, ~] = homerun(temp_state, temp_system, temp_experiment, temp_options, []);

        temp_eigenvalues(:, n) = eig(temp_system.ham);
        temp_integrals(:, n) = mean(real(detected_signal.sf), 2);
        LF_asymmetries(n, tau_idx) = real(detected_signal.sf(1, end));

        if save_all_data
            if save_traces_to_disk
                output_filename = sprintf('trace_k%03d_n%05d.mat', n, n);
                output_path = fullfile(save_trace_dir, output_filename);

                % sf = detected_signal.sf;        % time-domain signal

                % output.detected_signal = detected_signal.sf;
                % output.options = options;
                
                % Save traces to disk
                save(output_path, "-struct", 'output', '-v7.3');

                % Clear from worker's RAM
                % sf = [];
                detected_signal = [];
            else
                signals{n} = detected_signal.sf;
                % temp_allsignals{n} = detected_signal;
            end
        end
        
        eigenvalues{n} = temp_eigenvalues;
        % spectra(n,:,:) = temp_integrals;
    end

    det_op = 1;
    weights = sin(thetas);
    weights = weights / sum(weights);
    
    % Nt already set from experiment for this tau
    powder_signal_tau = zeros(Nt, 1);   % column vector
    
    % store time axis for this tau
    time_axes{tau_idx} = (0:Nt-1).' * experiment.dt;   % Nt x 1
    
    for n = 1:Norient
        % signals{n} is ndet x Nt, take row det_op and force column
        trace_n = signals{n}(det_op, :).';   % now Nt x 1 (transpose)
        powder_signal_tau = powder_signal_tau + weights(n) * trace_n;
    end
    
    powder_signals{tau_idx} = powder_signal_tau;   % Nt x 1
end

toc;

%% Plot sin(theta) weighted powder timetrace

fig = figure('NumberTitle','off', 'Name', 'Powder Averaged Time traces');
hold on
for tau_idx = 1:length(tau_array)
    plot(time_axes{tau_idx}(:), real(powder_signals{tau_idx}(:)));
end
hold off
xlabel('time (ns)')         % adjust units as needed
ylabel('P_z')
legend(arrayfun(@(t) sprintf('\\tau=%g ns', t), tau_array, 'UniformOutput', false));

% tau_idx = 3;

% plot(time_axes{tau_idx}, real(powder_signals{tau_idx}));

%% Plot LF asymmetries as a function of tau

fig = figure('NumberTitle','off', 'Name', 'LF Asymmetry');
hold on
scatter(tau_array, LF_asymmetries(1, :))
% legend_strings = arrayfun(@(x) sprintf('\tau = %.0f ns', x), tau_array, 'UniformOutput', false);
% legend(legend_strings)
hold off

%% Integrate over thetas
weights = sin(thetas);
det_op = 1;
powder_spectrum = zeros(Nfields, 1);
for ii = 1:Nfields
    powder_spectrum(ii) = sum(weights * spectra(:, det_op, ii))/sum(weights);
end

fig = figure('NumberTitle','off','Name','MQ Powder Spectrum');
hold on
plot(magnetic_fields, powder_spectrum)
xlabel('B / T')
ylabel('P_z')

% save('Data/num_ALC_simulation_SrTiO3_powdedwdwdwdr.mat', 'magnetic_fields', "powder_spectrum")

%% Plot time evolution of signal
% [experiment, options] = triple(sequence, options);  % build experiment to get experiment.tp
det_op = 1;

signal = zeros(size(signals{det_op}));
Nt = size(signals{det_op}, 2);
time = (0:Nt-1) * experiment.dt;
trace_idx = 1;

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

    % save('Data/MQ_signal_time_evolution_on_resonance.mat', 'trace', 'time_ds')
end

%% Plot MQ Spectra for different thetas
det_op = 1;

legendStrings = arrayfun(@(x) sprintf('\\theta = %.1f°', x), rad2deg(thetas), 'UniformOutput', false);

fig = figure('NumberTitle','off','Name','MQ different theta Spectra A_iso5'); 
hold on
h = gobjects(Norient,1);
for ii = 1:Norient
    h(ii) = plot(magnetic_fields, squeeze(spectra(ii, det_op, :)));
end
legend(h, legendStrings)
hold off

xlabel('B / T')
ylabel('P_z')

% Use print with -painters (vector graphics) and -dpdf
% exportgraphics(fig, 'C:\Users\walk_d\GitHub\MQ_muSR\Figures\ALC_simulations\Sim_different_theta.pdf', 'ContentType','vector','BackgroundColor','none')
