clear system
clear options
clear sequence

% Zeeman (gamma / GHz/T)
ge = 28.02495;
gmu = -0.1355;

% Coupling constants (in GHz) and rotation angles (in degree)
A_iso = 0.5148;
D_parallel = 0.002;
D_perpen = -D_parallel/2;

% thetas = deg2rad([1, 5, 20, 45, 70, 85, 89]);
% phis = deg2rad([0, 30, 45, 60, 90]);

thetas = deg2rad([5, 45]);
phis = deg2rad([45]);

% Range of B0
magnetic_fields = linspace(0, 6, 200);

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

%% Simulating for the first magnetic field in magnetid_fields
% nu_muon = gmu*magnetic_fields(400);
% nu_electron = ge*magnetic_fields(400);
% 
% % the system
% system.sqn=[0.5 0.5];       % spin quantum numbers
% system.interactions={1,0, 'z','e', nu_electron;       % Zeeman term of electron in simulation frame
%                      2,0, 'z','e', nu_muon;           % Zeeman of muon
%                      1,2, 'x','x', A_iso;             % isotropic HF
%                      1,2, 'y','y', A_iso;             % isotropic HF
%                      1,2, 'z','z', A_iso;             % isotropic HF
%                      1,2, 'x','x', D_perpen;          % axially symmetric HF
%                      1,2, 'y','y', D_perpen;          % axially symmetric HF
%                      1,2, 'z','z', D_parallel};       % axially symmetric HF
% 
% system.init_state='ez'; % LF mode (muon in the Iz eigenstate)
% system.eq = 'init';  % equilibrium state is the same as initial state
% 
% % the options
% options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
% options.down_conversion=0;  % downconversion of signal (1) or not (0)
% options.det_op={'ez', 'ex'};
% 
% options.labframe = 1;     % lab frame simulation is on
% options.awg.s_rate = 500;  % can be switched on to improve accuracy
% 
% sequence.tp=500.0;     % vector with event lengths in ns
% % sequence.nu1=[0];     % Pulse amplitude (no pulse when 0)
% % sequence.frq=[0, 0];
% sequence.detection=ones(1,length(sequence.tp)); % detection always on
% 
% [experiment,options] = triple(sequence, options); % builds experiment (pulses)
% 
% [system, state, options]=setup(system,options); % build system 
% 
% [state, detected_signal, experiment] = homerun(state, system, experiment, options, []);



%%
% the system
system.sqn=[0.5 0.5];       % spin quantum numbers

system.interactions = {};
                     
system.init_state='ez'; % LF mode (muon in the Iz eigenstate)
system.eq = 'init';  % equilibrium state is the same as initial state

% Set options
options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
options.det_op={'ez', 'ex'};
options.labframe = 1;     % lab frame simulation is on
options.awg.s_rate = 500;  % can be switched on to improve accuracy (gives sampling rate of simulation in GHz)

sequence.tp=500.0;     % vector with event lengths in ns
sequence.detection=ones(1,length(sequence.tp)); % detection always on

% Preallocate: one cell per orientation. Each cell holds a matrix: (#det_ops) x Nfields
Nfields = numel(magnetic_fields);
Norient = numel(T_labs);

signals = cell(1, Norient);
allsignals = cell(1, Norient);     % only if you need full detected_signal per field
eigenvalues = cell(1, Norient);

parfor k = 1:Norient
    temp_signal = zeros(length(options.det_op), Nfields);
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

        % Dipolar Hamiltonian in lab frame
        H_dip = T_lab(1,1)*Sx*Ix + T_lab(1,2)*Sx*Iy + T_lab(1,3)*Sx*Iz + ...
                T_lab(2,1)*Sy*Ix + T_lab(2,2)*Sy*Iy + T_lab(2,3)*Sy*Iz + ...
                T_lab(3,1)*Sz*Ix + T_lab(3,2)*Sz*Iy + T_lab(3,3)*Sz*Iz;

        temp_system.ham = temp_system.ham + H_dip;

        e = eig(temp_system.ham);
        temp_eigenvalues(:, n) = real(e);

        [temp_state, detected_signal, ~] = homerun(temp_state, temp_system, temp_experiment, temp_options, []);
        temp_signal(:, n) = detected_signal.sf;
        allsignals{k}{n} = detected_signal;
    end
    signals{k} = temp_signal;
    eigenvalues{k} = temp_eigenvalues;
end


toc;

[experiment, options] = triple(sequence, options);  % build experiment to get experiment variable

signal = zeros(size(signals{1}));
time = experiment.tp;

% Integrate
integrals = zeros(length(options.det_op), length(magnetic_fields));

% Filter is a lowpass butterworth filter of 3rd order
for n=1:length(magnetic_fields)
    signal = signal + signals{n};

    integrals(:, n) = sum(real(signals{n}), 2)/length(signal);
end

%%
% figure(7); clf; hold on
% plot(magnetic_fields, integrals(1, :))
% xlabel('Magnetic field B_0 / T')
% ylabel('Muon Polarization')


%% 
% Each column corresponds to one field, each row to an eigenvalue trajectory
E_matrix = cell2mat(eigenvalues);
E_matrix = reshape(E_matrix, size(eigenvalues{1},1), []); 

upper_diff = E_matrix(4, :) - E_matrix(3, :);
lower_diff = E_matrix(2, :) - E_matrix(1, :);

% Plot eigenvalues vs. field
figure(10); hold on
plot(magnetic_fields, real(E_matrix.'), 'LineWidth', 1.2)
xlabel('Magnetic field B_0 / T')
ylabel('Energy (GHz)')
title('Hamiltonian eigenvalues vs. magnetic field')

figure(11); hold on
plot(magnetic_fields, upper_diff, magnetic_fields, lower_diff)
xlabel('Magnetic field B_0 / T')
ylabel('Energy Difference / GHz')
legend('Top two levels', 'Bottom two levels')

mask = magnetic_fields > 0.1;

masked_magnetic_fields = magnetic_fields(mask);
masked_upper_diff = upper_diff(mask);

[min_value, min_index] = min(masked_upper_diff);
disp(masked_magnetic_fields(min_index))