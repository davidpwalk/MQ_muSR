clear system
clear options
clear sequence

%-- Settings --%
save_all_data = false;

% Zeeman (gamma / GHz/T)
ge = 28.02495;
gmu = -0.1355;

% Coupling constants (in GHz) and rotation angles (in degree)
A_iso = 0.5148;
D_parallel = 0.002;
D_perpen = -D_parallel/2;

thetas = deg2rad(linspace(0, 90, 200));
% thetas = deg2rad([1, 5, 20, 45, 70, 85, 89]);
phis = deg2rad([0]); % Phi has no impact on the spectra

% Range of B0
magnetic_fields = linspace(1.82, 1.98, 400);

% the system
system.sqn=[0.5 0.5];       % spin quantum numbers

system.interactions = {};
                     
system.init_state='ez'; % LF mode (muon in the Iz eigenstate)
system.eq = 'init';  % equilibrium state is the same as initial state

% Add relaxation
system.T1 = 2200;  % 2.2 us
system.T2 = 1e10;  % 10 s 

% Set options
options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex'};
options.labframe = 1;       % lab frame simulation is on
options.awg.s_rate = 12;   % gives sampling rate of simulation in GHz

sequence.tp=2000.0;     % vector with event lengths in ns
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

eigenvalues = cell(1, Norient);
spectra = zeros(Norient, numel(options.det_op), Nfields);

if save_all_data
    signals = cell(1, Norient);
    allsignals = cell(Norient, Nfields);     % only if you need full detected_signal per field
end

tic;

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

        % disp(H_dip)
        % disp('before')
        % disp(temp_system.ham)

        temp_system.ham = temp_system.ham + 2*pi*(H_dip + H_iso);
        
        % disp('By hand')
        % disp(temp_system.ham)

        temp_eigenvalues(:, n) = eig(temp_system.ham);

        [temp_state, detected_signal, ~] = homerun(temp_state, temp_system, temp_experiment, temp_options, []);

        if save_all_data
            temp_signal(:, :, n) = detected_signal.sf;
            temp_allsignals{n} = detected_signal;
        end

        temp_integrals(:,n) = mean(real(detected_signal.sf), 2);
    end
    
    eigenvalues{k} = temp_eigenvalues;
    spectra(k,:,:) = temp_integrals;

    if save_all_data
        signals{k} = temp_signal;
        allsignals(k,:) = temp_allsignals; 
    end
end

toc;

%%
% [experiment, options] = triple(sequence, options);  % build experiment to get experiment.tp
% 
% signal = zeros(size(signals{1}));
% Nt = size(signals{1}, 2);
% time = (0:Nt-1) * experiment.dt;
% 
% % Integrate
% integrals = zeros(Norient, numel(options.det_op), numel(magnetic_fields));
% 
% for k = 1:Norient
%     for n = 1:Nfields
%         % integrate over time (dimension 2)
%         integrals(k, 1, n) = mean(real(allsignals{k,n}.sf(1,:)));
%         integrals(k, 2, n) = mean(real(allsignals{k,n}.sf(2,:)));
%         % squeeze(sum(real(signals{k}), 2)) / Nt;
%     end
% end

%% Plot time evolution of signal
% stride = 100;   % downsample
% t_idx = 1:stride:length(time);
% t_idx = 1:stride;
% trace = squeeze(signals{1}(1, t_idx, end));
% time_ds = time(t_idx);
% 
% figure(343); clf; hold on
% plot(time_ds, real(trace))
% xlabel('Time / ns')
% ylabel('Polarization')

%% 
E_matrix = cell2mat(eigenvalues);

% reshape to 3-D: (4 eigenvals) × (Nfields) × (Norient)
E_matrix_3D = reshape(E_matrix, size(eigenvalues{1},1), Nfields, Norient);

k = 1;
figure(10); clf; hold on
plot(magnetic_fields, real(squeeze(E_matrix_3D(:,:,k))').', 'LineWidth', 1.2)  % -> Nfields x 4
xlabel('Magnetic field B_0 / T')
ylabel('Energy (GHz)')
title(sprintf('Eigenvalues vs B (orientation %d)', k))

k = 1;
figure(10); clf; hold on
plot(magnetic_fields, real(squeeze(E_matrix_3D(:,:,k))').', 'LineWidth', 0.8)
xlabel('Magnetic field B_0 / T')
ylabel('Energy (GHz)')

k = 1;
upper_diff = real(E_matrix_3D(4,:,k) - E_matrix_3D(3,:,k));   % 1 × Nfields
lower_diff = real(E_matrix_3D(2,:,k) - E_matrix_3D(1,:,k));   % 1 × Nfields

figure(11); clf; hold on
plot(magnetic_fields, upper_diff, magnetic_fields, lower_diff)
xlabel('Magnetic field B_0 / T')
ylabel('Energy Difference / GHz')
legend('Top two levels', 'Bottom two levels')

mask = magnetic_fields > 0.1;

masked_magnetic_fields = magnetic_fields(mask);
masked_upper_diff = upper_diff(mask);

[min_value, min_index] = min(masked_upper_diff);
level_crossing = masked_magnetic_fields(min_index);
disp(level_crossing)

%% get the traces in spectral domain (only works if save_all_data = true)

if save_all_data
    tdy = zeros(length(allsignals{1}.signals.t{end}),Nfields);
    tdx = allsignals{1}.signals.t{1};
    
    for ii=1:Nfields
    %     tdy(:,ii) = real(signalc{ii}(1,:));
    %     tdy(:,ii) = allsignal{ii}.signals.dc{end}(3,:);
        % Ix
    %     tdy(:,ii) = allsignal{ii}.signals.sf{end}(4,:);
    %     tdy(:,ii) = tdy(:,ii) - mean(tdy(:,ii));
        
        % Iz
        tdy(:,ii) = allsignals{1,ii}.signals.sf{1}(1,:);
    
    
        % Ix
        tdy(:,ii) = allsignals{1,ii}.signals.sf{1}(2,:);
        
        
    end
    
    
    % get spectra
    [fdy, fdx] = getmultspec(tdy,tdx,2,1);
    
    fdx = fdx;
    p_idx = fdx < 100 & fdx > -100;
    p_idx = 1:length(fdx);
    
    i2Dcut_nox(fdx(p_idx),magnetic_fields,abs(fdy(p_idx,:)),23);
    %xlim(0,100);
    title('Hand Hamiltonian')
    
    figure(200);
    clf();
    contour(fdx(p_idx),magnetic_fields,abs(fdy(p_idx,:))',22)
    % ylim([1370, 1430])
    % xlim([-80, 80])
    xlabel('f [GHz]')
end
%% Plot ALC Spectra for different thetas
det_op = 1;
peak_positions = [];

% Create figure
fig = figure('NumberTitle','off','Name','ALC Spectra'); 
hold on
for ii = 1:Norient
    [~, min_index] = min(spectra(ii, det_op, :));
    peak_positions(ii) = magnetic_fields(min_index);
    plot(magnetic_fields, squeeze(spectra(ii, det_op, :)))
end
peak_positions = [rad2deg(thetas(:)), peak_positions(:)];
peak_positions(peak_positions(:, 2) < 1.85, :) = []; % Drop thetas where there is no peak at all
hold off
xlabel('B / T')
ylabel('P_z')
% xlim([1.8675, 1.9325])

% save('Data/num_ALC_simulation_thetas', 'magnetic_fields', "spectra")

legendStrings = arrayfun(@(x) sprintf('\\theta = %.1f°', x), rad2deg(thetas), 'UniformOutput', false);
legend(legendStrings, 'Location', 'best')

% Use print with -painters (vector graphics) and -dpdf
% exportgraphics(fig, 'C:\Users\walk_d\GitHub\MQ_muSR\Figures\ALC_simulations\Sim_different_theta.pdf', 'ContentType','vector','BackgroundColor','none')

% fig = figure('NumberTitle','off','Name','Peak Positions vs theta');
% hold on;
% plot(peak_positions(:, 1), peak_positions(:, 2))
% hold off;
%% Integrate over thetas

weights = sin(thetas);
powder_spectrum = zeros(Nfields, 1);
% TODO: check if the weighted sum is actually calculated correctly
for ii = 1:Nfields
    powder_spectrum(ii) = sum(weights * spectra(:, det_op, ii))/sum(weights);
end

fig = figure('NumberTitle','off','Name','ALC Powder Spectrum');
hold on
plot(magnetic_fields, powder_spectrum)
xlabel('B / T')
ylabel('P_z')

save('Data/num_ALC_simulation_powder', 'magnetic_fields', "powder_spectrum")

%% Compare peak positions to analytical counterpart

s = load("ana_peak_positions.mat");
ana_peak_positions = s.peak_positions;

peak_positions_diff = ana_peak_positions(:, 2) - peak_positions(:, 2);

fig = figure('NumberTitle','off','Name','Peak Pos difference');
hold on
% plot(peak_positions(:, 1), peak_positions_diff)
hold off

save('peak_positions_0_5awg.mat','peak_positions', 'peak_positions_diff')

filename = sprintf('results_all_thetas_awg%.2f_tp%.2f.mat', options.awg.s_rate, sequence.tp);
% save(filename, 'options', 'system', 'sequence', 'peak_positions', 'peak_positions_diff')