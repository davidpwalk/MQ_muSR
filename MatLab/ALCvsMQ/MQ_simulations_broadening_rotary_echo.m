clear system
clear options
clear sequence

%-- Settings --%
rot_echo = 0;
pulse_delay = 1000;  % in ns

% if save_traces_to_disk && ~exist("save_trace_dir", 'dir')
%     mkdir(save_trace_dir)
% end

% Sweep (default values used, while other value is sweeped)
nu1_default = 1;
T2_default = 2000;
tau_default = 200;

sweep_param = 'tau';
sweep_values = [200, 400, 800];

% Zeeman (gamma / GHz/T)
ge = -28.02495;
gmu = 0.1355;

% Coupling constants (in GHz) and rotation angles (in degree)
A_iso = 0.0014;
D_parallel = 0.0155;
% D_parallel = 0;
D_perpen = -D_parallel/2;

thetas = deg2rad(linspace(0, 90, 50));
% thetas = deg2rad([1, 5, 20, 45, 70, 85, 89]);
% thetas = deg2rad([0]);
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

% Set options
options.relaxation=1;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex', 'ze'};
options.labframe = 1;       % lab frame simulation is on
options.awg.s_rate = 20;   % gives sampling rate of simulation in GHz

nu_uw = abs(nu_electron);
% nu_uw = 3000

if rot_echo
    phase_flip = pi;
    sequence.phase = [0 ,0, phase_flip];  % last element will be overwritten in loop
    % sequence.nu1  = [0, 0, 0];
    sequence.frq  = [0, nu_uw, nu_uw];
    sequence.t_rise = [0, 0, 0];               % no chirp
else
    sequence.frq = [0, nu_uw];
    sequence.t_rise = [0, 0];
    sequence.phase = [0, 0];
end
sequence.type={'rectangular'};

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

spectra_taus = zeros(Norient, length(sweep_values));
LF_asymmetries = zeros(Norient, length(sweep_values));

% parpool('Threads');  % Start ThreadPool, threads share memory
for sweep_idx = 1:length(sweep_values)
    % fprintf('\n\nsweep_idx %d\n', sweep_idx)
    % disp('====================================================')
    
    if strcmp(sweep_param, 'T2')
        T2 = sweep_values(sweep_idx);
        nu1 = nu1_default;
        tau = tau_default;
        unit = 'ns';  % unit displayed in plot legend
    
    elseif strcmp(sweep_param, 'nu1')
        nu1 = sweep_values(sweep_idx);
        T2 = T2_default;
        tau = tau_default;
        unit = 'MHz';  % unit displayed in plot legend
    elseif strcmp(sweep_param, 'tau')
        tau = sweep_values(sweep_idx);
        nu1 = nu1_default;
        T2 = T2_default;
        unit = 'ns';
    end

    system.T2 = [0, 1e9, T2,  T2;
                 0,   0, T2,  T2;            
                 0,   0,  0, 1e9;
                 0,   0,  0,   0];

    if rot_echo
        sequence.phase(3) = phase_flip+angle(exp(1j*2*pi*nu_uw*(tau)));  % account for acquired phase
        sequence.tp = [pulse_delay, tau, tau];  % pulse delay with two same pulse blocks afterwards
        sequence.nu1  = [0, nu1, nu1];
    else
        sequence.tp = [pulse_delay, 2*tau];  % pulse_delay + pulse with length 2*tau
        sequence.nu1 = [0, nu1];
    end

    % Testing
    sequence.nu1 = [0, 1];

    sequence.detection = ones(1,length(sequence.tp));
    
    [experiment, options] = triple(sequence, options);  % build experiment to get experiment.tp
    
    Nt = sum(experiment.tp)/sum(experiment.dt) + 1;

    % if save_all_data
    %     signals = cell(1, Norient);   % will hold ndet x Nt arrays
    % end

    parfor n = 1:Norient
        T_lab = T_labs{n};

        % temp_integrals = zeros(numel(options.det_op), Nfields);
    
        % if save_all_data
        %     temp_signal = zeros(length(options.det_op), Nt, Norient);
        %     % temp_allsignals = cell(1);
        % end
    
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

        temp_integrals(:, n) = mean(real(detected_signal.sf), 2);
        LF_asymmetries(n, sweep_idx) = real(detected_signal.sf(1, end));

        signals{n} = detected_signal.sf;
    end

    det_op = 1;
    weights = sin(thetas);
    weights = weights / sum(weights);
    
    % Nt already set from experiment for this tau
    powder_signal_tau = zeros(Nt, 1);   % column vector
    
    % store time axis for this tau
    time_axes{sweep_idx} = (0:Nt-1).' * experiment.dt;   % Nt x 1
    
    for n = 1:Norient
        % signals{n} is ndet x Nt, take row det_op and force column
        trace_n = signals{n}(det_op, :).';   % now Nt x 1 (transpose)
        powder_signal_tau = powder_signal_tau + weights(n) * trace_n;
    end
    
    powder_signals{sweep_idx} = powder_signal_tau;   % Nt x 1
end

toc;

%% Plot sin(theta) weighted powder timetrace

fig = figure('NumberTitle','off', 'Name', sprintf('Powder Averaged Time traces, delay=%d ns, rot_echo=%d', pulse_delay, rot_echo));
hold on
for sweep_idx = 1:length(sweep_values)
    plot(time_axes{sweep_idx}(:), real(powder_signals{sweep_idx}(:)));
end
hold off
xlabel('time (ns)')         % adjust units as needed
ylabel('P_z')
legend(arrayfun(@(t) sprintf('%s=%g %s', sweep_param, t, unit), sweep_values, 'UniformOutput', false));

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
