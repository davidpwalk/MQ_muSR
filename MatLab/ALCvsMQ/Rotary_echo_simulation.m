clear system
clear options
clear sequence

%-- Settings --%
save_all_data = false;

% Zeeman (gamma / GHz/T)
ge = -28.02495;
gmu = 0.1355;

% Coupling constants (in GHz) and rotation angles (in degree)
A_iso = 0.003;
D_parallel = 0.0155;
D_perpen = -D_parallel/2;

% thetas = deg2rad(linspace(0, 90, 20));
% thetas = deg2rad([1, 5, 20, 45, 70, 85, 89]);
thetas = deg2rad([45]);
phis = deg2rad([0]); % Phi has no impact on the spectra

magnetic_fields = [0.0957];

B0 = 0.0957;

% B1 inhomogenity
B1_center = 0.0822;
B1_sigma = 0;
B1_start = B1_center - 3*B1_sigma;
B1_end = B1_center + 3*B1_sigma;


% nu1 inhomogenity
nu1_center = 10;
nu1_sigma = 0.5;
nu1_start = nu1_center - 3*nu1_sigma;
nu1_end = nu1_center + 3*nu1_sigma;
nu1_step_size = 0.3;
nu1_vec = nu1_start:nu1_step_size:nu1_end;
nu1_probs = exp(-((nu1_center-nu1_vec)/nu1_sigma).^2);
nu1_probs = nu1_probs/trapz(nu1_probs);


% the system
system.sqn=[0.5 0.5];       % spin quantum numbers

system.interactions = {};
                     
system.init_state='ez'; % LF mode (muon in the Iz eigenstate)
system.eq = 'init';  % equilibrium state is the same as initial state

% Relaxtion matrix (e- transitions 2 us, rest 1 s)
system.T1 = [0, 1e9,  1e9,  1e9;
             0,   0,  1e9,  1e9;            
             0,   0,    0,  1e9;
             0,   0,    0,    0];

system.T2 = [0, 1e9, 2000,  1e9;
             0,   0,  1e9, 2000;            
             0,   0,    0,  1e9;
             0,   0,    0,    0];

% Set options
options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex', 'ze'};
options.labframe = 1;       % lab frame simulation is on
options.awg.s_rate = 6;   % gives sampling rate of simulation in GHz

nu_muon = gmu*B0;
nu_electron = ge*B0;

nu_uw = abs(nu_electron);
% nu_uw = 3000

% Set sequence
tau = 8000;  % ns, length of each pulse block

sequence.tp   = [tau, tau];             % two equal pulse blocks
sequence.nu1  = [1.0, 1.0];
sequence.phase = [0, pi];               % phase shift for second part of pulse
sequence.frq  = [nu_uw, nu_uw];
sequence.t_rise = [0, 0];               % no chirp
sequence.detection = ones(1,length(sequence.tp));

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
Nt = sum(experiment.tp)/sum(experiment.dt) + 1;

% Array for different nu2
inhom_array = cell(1, length(nu1_vec));

signals = cell(Nfields, Norient);

tic;

% parpool('Threads');  % Start ThreadPool, threads share memory
for k = 1:Nfields
    nu_muon = gmu * magnetic_fields(k);
    nu_electron = ge * magnetic_fields(k);

    % These are overwritten in the parfor loop anyway
    system.interactions = { ...
            1,0,'z','e', nu_electron; ...
            2,0,'z','e', nu_muon; ...
        };

    T_lab = T_labs{k};

    temp_signal = zeros(length(options.det_op), Nt, Nfields);
    temp_allsignals = cell(1, Nfields);

    for n = 1:Norient
        temp_T_lab = T_labs{n};

            
        % Isotropic Dipolar Hamiltonian
        H_iso = Sx*A_iso*Ix + Sy*A_iso*Iy + Sz*A_iso*Iz;
        
        % Dipolar Hamiltonian in lab frame
        H_dip = Sx*T_lab(1,1)*Ix + Sx*T_lab(1,2)*Iy + Sx*T_lab(1,3)*Iz + ...
            Sy*T_lab(2,1)*Ix + Sy*T_lab(2,2)*Iy + Sy*T_lab(2,3)*Iz + ...
            Sz*T_lab(3,1)*Ix + Sz*T_lab(3,2)*Iy + Sz*T_lab(3,3)*Iz;
       
        [experiment, options] = triple(sequence, options);
        [system, ~, ~] = setup(system, options);

        system.ham = system.ham + 2*pi*(H_dip + H_iso);

        for i = 1:length(nu1_vec)
            temp_system = system;
            temp_system.interactions{1, end} = nu1_vec(i);

            [temp_system, temp_state, temp_options] = setup(temp_system, options);
    
            [temp_state, detected_signal, ~] = homerun(temp_state, temp_system, experiment, temp_options, []);
    
            inhom_array{i} = detected_signal.sf;
        end
        signal = zeros(size(inhom_array{1}));

        for i = 1:length(nu1_vec)
            signal = signal + inhom_array{i}*nu1_probs(i);
        end

        % Store signal for one field and orientation but with added
        % inhomogenity
        signals{k, n} = signal;
    end
end

toc;

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
det_op = 1;
orientation = 1;
magnetic_field = 1;

time = (0:Nt-1) * experiment.dt;

stride = 1;   % downsample
t_idx = 1:stride:length(time);
trace = signals{magnetic_field, orientation}(det_op, :);
time_ds = time(t_idx);

fig = figure('NumberTitle','off','Name','Time-Domain Spectrum Pz');
hold on
plot(time_ds, real(trace))
xlabel('Time / ns')
ylabel('Polarization')

% save('Data/MQ_signal_timdwdwdwde_evolution_on_resonance.mat', 'trace', 'time_ds')

%% Plot MQ Spectra for different thetas
det_op = 1;

legendStrings = arrayfun(@(x) sprintf('\\theta = %.1f°', x), rad2deg(thetas), 'UniformOutput', false);

fig = figure('NumberTitle','off','Name','MQ different theta Spectra'); 
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
