clear system
clear options
clear sequence

% Zeeman (gamma / GHz/T)
ge = 28.02495;
gmu = -0.1355;

% Coupling constants (in GHz) and rotation angles (in degree)
A_iso = 0.5148;
A_iso = 0.0;
A_iso = 0.02;
D_parallel = 0.002;
% D_parallel = 0.020;
D_parallel = 0.0;
D_perpen = -D_parallel/2;

% thetas = deg2rad([1, 5, 20, 45, 70, 85, 89]);
% phis = deg2rad([0, 30, 45, 60, 90]);

thetas = deg2rad([45]);
phis = deg2rad([45]);

% Range of B0
magnetic_fields = linspace(0, 6, 200);
% magnetic_fields = linspace(1.8, 2.0, 200);
% magnetic_fields = linspace(1, 3, 200);
% magnetic_fields = [4, 5 , 6];
magnetic_fields = linspace(0, 6, 20);

Sx = sop([1/2 1/2],'xe');
Sy = sop([1/2 1/2],'ye');
Sz = sop([1/2 1/2],'ze');

Ix = sop([1/2 1/2],'ex');
Iy = sop([1/2 1/2],'ey');
Iz = sop([1/2 1/2],'ez');

% Rotation matrices
Rz = @(theta) [cos(theta) -sin(theta) 0;
               sin(theta)  cos(theta) 0;
               0           0          1];

Ry = @(phi) [cos(phi)  0 sin(phi);
             0         1        0;
             -sin(phi) 0 cos(phi)];

T_principal = diag([D_perpen, D_perpen, D_parallel]);

% Number of combinations
N = numel(thetas) * numel(phis);

% Preallocate cell array
T_labs = cell(1, N);

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

%%
% the system
system.sqn=[0.5 0.5];       % spin quantum numbers

system.interactions = {};
                     
system.init_state='ex'; % LF mode (muon in the Iz eigenstate)
system.eq = 'init';  % equilibrium state is the same as initial state

% Set options
options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex'};
options.labframe = 1;     % lab frame simulation is on
options.awg.s_rate = 500;  % can be switched on to improve accuracy (gives sampling rate of simulation in GHz)

sequence.tp=500.0;     % vector with event lengths in ns
sequence.detection=ones(1,length(sequence.tp)); % detection always on

[experiment, options] = triple(sequence, options);  % build experiment to get experiment.tp

% Preallocate: one cell per orientation. Each cell holds a matrix: (#det_ops) x Nfields
Nfields = numel(magnetic_fields);
Norient = numel(T_labs);
Nt = numel(experiment.tp);

signals = cell(1, Norient);
allsignals = cell(Norient, Nfields);     % only if you need full detected_signal per field
eigenvalues = cell(1, Norient);

for k = 1:Norient
    temp_eigenvalues = zeros(4, Nfields);

    T_lab = T_labs{k};

    for n = 1:Nfields
        temp_nu_muon = gmu * magnetic_fields(n);
        temp_nu_electron = ge * magnetic_fields(n);

        temp_system = system;  % copy
        temp_system.interactions = { ...
            1,0,'z','e', temp_nu_electron; ...
            2,0,'z','e', temp_nu_muon; ...
            1,2,'x','x', A_iso; ...
            1,2,'y','y', A_iso; ...
            1,2,'z','z', A_iso ...
        };

        [temp_experiment, temp_options] = triple(sequence, options);
        [temp_system, temp_state, temp_options] = setup(temp_system, temp_options);

        % Isotropic dipolar Hamiltonian
        A_iso_matrix = 2*pi*diag([A_iso, A_iso, A_iso]);
        H_iso = Sx*A_iso_matrix(1, 1)*Ix + Sy*A_iso_matrix(2, 2)*Iy + Sz*A_iso_matrix(3, 3)*Iz;

        % Dipolar Hamiltonian in lab frame
        H_dip = Sx*T_lab(1,1)*Ix + Sx*T_lab(1,2)*Iy + Sx*T_lab(1,3)*Iz + ...
                Sy*T_lab(2,1)*Ix + Sy*T_lab(2,2)*Iy + Sy*T_lab(2,3)*Iz + ...
                Sz*T_lab(3,1)*Ix + Sz*T_lab(3,2)*Iy + Sz*T_lab(3,3)*Iz;

        % disp(H_dip)
        % 
        % disp('before')
        % disp(temp_system.ham)

        temp_system.ham = temp_system.ham + 2*pi*(H_dip);
        
        % disp('By hand')
        % disp(temp_system.ham)

        temp_eigenvalues(:, n) = eig(temp_system.ham);

        [temp_state, detected_signal, ~] = homerun(temp_state, temp_system, temp_experiment, temp_options, []);

        Nt = size(detected_signal.sf, 2);
        temp_signal = zeros(length(options.det_op), Nt, Nfields);
        temp_signal(:, :, n) = detected_signal.sf;
        allsignals{k,n} = detected_signal;
    end
    signals{k} = temp_signal;
    eigenvalues{k} = temp_eigenvalues;
end

toc;

%%
[experiment, options] = triple(sequence, options);  % build experiment to get experiment.tp

signal = zeros(size(signals{1}));
Nt = size(signals{1}, 2);
time = (0:Nt-1) * experiment.dt;

% Integrate
integrals = zeros(Norient, numel(options.det_op), numel(magnetic_fields));

for k = 1:Norient
    for n = 1:Nfields
        % integrate over time (dimension 2)
        integrals(k, 1, n) = mean(real(allsignals{k,n}.sf(1,:)));
        integrals(k, 2, n) = mean(real(allsignals{k,n}.sf(2,:)));
        % squeeze(sum(real(signals{k}), 2)) / Nt;
    end
end

%%
stride = 100;   % downsample
t_idx = 1:stride:length(time);
t_idx = 1:stride;
trace = squeeze(signals{1}(1, t_idx, end));
time_ds = time(t_idx);

figure(343); clf; hold on
plot(time_ds, real(trace))
xlabel('Time / ns')
ylabel('Polarization')

%%
figure(7); clf; hold on
plot(magnetic_fields, squeeze(integrals(1, 1, :)))
xlabel('Magnetic field B_0 / T')
ylabel('<I_z>')


%% 
% Each column corresponds to one field, each row to an eigenvalue trajectory
E_matrix = cell2mat(eigenvalues);
E_matrix = reshape(E_matrix, size(eigenvalues{1},1), []); 

upper_diff = E_matrix(4, :) - E_matrix(3, :);
lower_diff = E_matrix(2, :) - E_matrix(1, :);

% Plot eigenvalues vs. field
figure(10); clf; hold on
plot(magnetic_fields, real(E_matrix.'), 'LineWidth', 1.2)
xlabel('Magnetic field B_0 / T')
ylabel('Energy (GHz)')
title('Hamiltonian eigenvalues vs. magnetic field')

figure(11); clf; hold on
% plot(magnetic_fields, upper_diff, magnetic_fields, lower_diff)
plot(magnetic_fields, upper_diff)
xlabel('Magnetic field B_0 / T')
ylabel('Energy Difference / GHz')
legend('Top two levels', 'Bottom two levels')

mask = magnetic_fields > 0.1;

masked_magnetic_fields = magnetic_fields(mask);
masked_upper_diff = upper_diff(mask);

[min_value, min_index] = min(masked_upper_diff);
disp(masked_magnetic_fields(min_index))



%% get the traces in spectral domain


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

i2Dcut_nox(fdx(p_idx),magnetic_fields,abs(fdy(p_idx,:)),25);
%xlim(0,100);
title('Spidyan Hamiltonian')

figure(200);
clf();
contour(fdx(p_idx),magnetic_fields,abs(fdy(p_idx,:))',24)
% ylim([1370, 1430])
% xlim([-80, 80])
xlabel('f [GHz]')

