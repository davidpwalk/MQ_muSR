clear system
clear options
clear sequence

%-- Settings --%
save_all_data = true;
save_traces_to_disk = false;  % only applies if save_all_data=true
save_trace_dir = 'Data/traces/test';

if save_traces_to_disk && ~exist("save_trace_dir", 'dir')
    mkdir(save_trace_dir)
end

% Sweep settings
nu1_default = 1;  % value used if other param is sweeped
T2_default = 2000;

sweep_param = 'nu1';
% sweep_values = [5, 10, 20, 50, 100, 1000, 2000];
sweep_values = [1, 2, 10, 50, 100];

% Zeeman (gamma / GHz/T)
ge = -28.02495;
gmu = 0.1355;

% Coupling constants (in GHz) and rotation angles (in degree)
A_iso = 0.0014;
D_parallel = 0.0155;
D_perpen = -D_parallel/2;

% thetas = deg2rad(linspace(0, 90, 50));
% thetas = deg2rad([1, 5, 20, 45, 70, 85, 89]);
thetas = deg2rad([45]);
phis = deg2rad([0]); % Phi has no impact on the spectra

% Range of B0
B0 = 0.0822;
B_start = 0.078;
B_end = 0.086;
dB = 0.0001;
magnetic_fields = B_start : dB : B_end;
% magnetic_fields = linspace(0, 0.16, 200);
magnetic_fields = [0.79, 0.08177];

length(magnetic_fields)

% the system
system.sqn=[0.5 0.5];       % spin quantum numbers

system.interactions = {};
                     
system.init_state='ez'; % LF mode (muon spin polarized)
% system.init_state='zz'; % both e- and muon in Jz eigenstate (only SzIz occupied)
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
options.awg.s_rate = 6;   % gives sampling rate of simulation in GHz

%

nu_muon = gmu*B0;
nu_electron = ge*B0;

nu_uw = abs(nu_electron);
% nu_uw = 3000

% Set sequence
sequence.tp=16000.0;                             % vector with event lengths in ns
sequence.frq=nu_uw;                              % frequency of pulse
sequence.t_rise=0;                               % rise time of chirp pulses
sequence.detection=ones(1,length(sequence.tp));  % detection always on

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

eigenvalues = cell(1, Norient);
spectra = zeros(Norient, numel(options.det_op), Nfields);

if save_all_data
    signals = cell(1, Norient);
    allsignals = cell(Norient, Nfields);     % only if you need full detected_signal per field
end

tic;

spectra_sweep = zeros(Norient, Nfields, length(sweep_values));

for sweep_idx = 1:length(sweep_values)

if strcmp(sweep_param, 'T2')
    T2 = sweep_values(sweep_idx);
    nu1 = nu1_default;
    unit = 'ns';  % unit displayed in plot legend

elseif strcmp(sweep_param, 'nu1')
    nu1 = sweep_values(sweep_idx);
    T2 = T2_default;
    unit = 'MHz';  % unit displayed in plot legend
end

sequence.nu1 = nu1;
system.T2 = [0, 1e9, T2,  T2;
             0,   0, T2,  T2;            
             0,   0,  0, 1e9;
             0,   0,  0,   0];

% parpool('Threads');  % Start ThreadPool, threads share memory
for k = 1:Norient
    T_lab = T_labs{k};
    temp_eigenvalues = zeros(4, Nfields);
    temp_integrals = zeros(numel(options.det_op), Nfields);

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

        temp_system.ham = temp_system.ham + 2*pi*(H_dip + H_iso);  % 2pi needed cos they are in angular frequency units in spidyan

        [temp_state, detected_signal, ~] = homerun(temp_state, temp_system, temp_experiment, temp_options, []);

        temp_eigenvalues(:, n) = eig(temp_system.ham);
        temp_integrals(:,n) = mean(real(detected_signal.sf), 2);

        if save_all_data
            if save_traces_to_disk
                output_filename = sprintf('trace_k%03d_n%05d.mat', k, n);
                output_path = fullfile(save_trace_dir, output_filename);

                % sf = detected_signal.sf;        % time-domain signal

                % output.detected_signal = detected_signal.sf;
                % output.options = options;
                
                % Save traces to disk
                save(output_path, "-struct", 'output', '-v7.3');

                % Clear from worker’s RAM
                % sf = [];
                detected_signal = [];
            else
                temp_signal(:, :, n) = detected_signal.sf;
                temp_allsignals{n} = detected_signal;
            end
        end
    end
    
    eigenvalues{k} = temp_eigenvalues;
    spectra(k,:,:) = temp_integrals;

    if save_all_data
        signals_sweep{sweep_idx} = zeros(numel(options.det_op), Nt, Nfields, Norient);
        for k = 1:Norient
            signals_sweep{sweep_idx}(:,:,:,k) = signals{k};
        end
    end
end

det_op = 1;
spectra_sweep(:, :, sweep_idx) = squeeze(spectra(:, det_op, :));

end

toc;

%% Show spectra for sweep
fig = figure('NumberTitle','off','Name', join('MQ spectra different ', sweep_param));
hold on
for ii = 1:length(sweep_values)
    plot(magnetic_fields, spectra_sweep(:, :, ii))
end
legend_strings = arrayfun(@(x) sprintf('%s = %.0f %s', sweep_param, x, unit), sweep_values, 'UniformOutput', false);
legend(legend_strings)
hold off

%% Integrate over thetas
weights = sin(thetas);
det_op = 1;
powder_spectra = cell(1, length(sweep_values));
norm_weights = weights(:) / sum(weights);

for sweep_idx = 1:length(sweep_values)
    powder_spectra{sweep_idx} = (norm_weights.' * spectra_sweep(:, :, sweep_idx)).';   % (Nfields x 1)
end

fig = figure('NumberTitle','off','Name','MQ Powder Spectrum');
hold on
for sweep_idx = 1:length(sweep_values)
    plot(magnetic_fields, powder_spectra{sweep_idx}, 'DisplayName', sprintf('T2_{elec} = %.0f ns', sweep_values(sweep_idx)))
end
legend show
xlabel('B / T')
ylabel('P_z')

% save('Data/num_MQ_powder_spectra_D15_5MHz_different_T2dwdadawd.mat', 'magnetic_fields', "powder_spectra", "T2_array")

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

%% Time-domain signal with magnetic-field slider

det_op = 1;
orientation = 1;
sweep_idx = 1;

time = (0:Nt-1) * experiment.dt;
trace_all = squeeze(signals{orientation}(det_op, :, :)); % Nt × Nfields

fig = figure('NumberTitle','off','Name','Time-Domain Signal');
ax = axes(fig); hold(ax,'on');

plt = plot(ax, time, real(trace_all(:,1)));
xlabel(ax,'Time / ns');
ylabel(ax,'Polarization');
title(ax, sprintf('B = %.5f T', magnetic_fields(1)));

sld = uicontrol(fig,'Style','slider', ...
    'Min',1,'Max',numel(magnetic_fields), ...
    'Value',1, ...
    'SliderStep',[1/(numel(magnetic_fields)-1), 5/(numel(magnetic_fields)-1)], ...
    'Units','normalized','Position',[0.15 0.02 0.7 0.04]);

sld.Callback = @(src,~) ...
    set(plt,'YData',real(trace_all(:,round(src.Value)))) | ...
    title(ax, sprintf('B = %.5f T', magnetic_fields(round(src.Value))));

%% Time-domain signal with sweep-value slider

det_op = 1;
orientation = 1;
field_idx = 1;        % <<< fixed magnetic-field index

time = (0:Nt-1) * experiment.dt;

% initial sweep
sweep_idx = 1;
trace_all = squeeze(signals_sweep{sweep_idx}{orientation}(det_op,:,field_idx));

fig = figure('NumberTitle','off','Name','Time-Domain Signal');
ax = axes(fig); hold(ax,'on');

plt = plot(ax, time, real(trace_all));
xlabel(ax,'Time / ns');
ylabel(ax,'Polarization');
title(ax, sprintf('%s = %.3g', sweep_param, sweep_values(sweep_idx)));

sld = uicontrol(fig,'Style','slider', ...
    'Min',1,'Max',numel(sweep_values), ...
    'Value',1, ...
    'SliderStep',[1/(numel(sweep_values)-1), 1/(numel(sweep_values)-1)], ...
    'Units','normalized','Position',[0.15 0.02 0.7 0.04]);

sld.Callback = @(src,~) ...
    set(plt,'YData', ...
        real(squeeze( ...
            signals_sweep{round(src.Value)}{orientation}(det_op,:,field_idx)))) | ...
    title(ax, sprintf('%s = %.3g', ...
        sweep_param, sweep_values(round(src.Value))));



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

% save('Data/SQ_pcolor_electrondwdawdw.mat', 'magnetic_fields', 'thetas', 'spectra')
