clear system
clear options
clear sequence



% Zeeman (gamma / GHz/T or MHz/mT)
ge = 28.02495;
% ge = 28.02495*0.98;
gmu = -0.1355;

% coupling
A_iso = 4.4956;

% magnetic field / T
B0 = 0.0822;

nu_muon = gmu*B0;
nu_electron = ge*B0;


% analytical solutions
nu_mid = (nu_electron + nu_muon)/2.0;
nu_diff = (nu_electron - nu_muon)/2.0;


alpha = 0.5*atan2(A_iso,nu_electron-nu_muon);
S =  sqrt((nu_diff)^2.0 + (A_iso/2)^2.0);
w12 = nu_mid + A_iso/2 - S;
w13 = nu_mid + A_iso/2 + S;
w24 = nu_mid - A_iso/2 + S;
w34 = nu_mid - A_iso/2 - S;

i12 = (1-sin(2*alpha))/4;
i13 = (1+sin(2*alpha))/4;
i24 = (1+sin(2*alpha))/4;
i34 = (1-sin(2*alpha))/4;


E1 = nu_mid + A_iso/2;
E2 = -A_iso/2 + S;
E3 = -A_iso/2 - S;
E4 = -nu_mid + A_iso/2;


% microwave frequency % Q?: change this?
nu_uw = abs(w34);

% vectorize (for fieldsweeps)

B_vec = [B0];


% add custom detection operators (x, y, z: component; a, b which spin)
d34 = sop([1/2 1/2],'bz');
d12 = sop([1/2 1/2],'az');
d13 = sop([1/2 1/2],'za');
d24 = sop([1/2 1/2],'zb');


% the zero quantum coherence, both as detection operator and for
% transformations
% I is the muon, S the e-
SxIy = sop([1/2 1/2],'xy');
SyIx = sop([1/2 1/2],'yx');

zq = 2*SxIy - 2*SyIx;

unitary = expm(-1i*alpha * zq);
unitary_t = expm(1i*alpha * zq);

Iz = sop([1/2 1/2],'ez'); % e for identiy, basically 1 tensordot I_z (expanding base)
Sz = sop([1/2 1/2],'ze');


Ix = sop([1/2 1/2],'ex');
Sx = sop([1/2 1/2],'xe');

% put detection operators from eigenframe to labframe
Iz_t = unitary_t * Iz * unitary;
d34_t = unitary_t * d34 * unitary;
d12_t = unitary_t * d12 * unitary;
d13_t = unitary_t * d13 * unitary;
d24_t = unitary_t * d24 * unitary;

% construct correct excitation operator with corresponding scalings
exc_op = ge * Sx + gmu * Ix; % nu1 now has units of MHz/mT, where mT are in the rotating frame
exc_op = exc_op/2; % now nu1 has units of mT in the laboratory frame (there is a factor of two in wave_propagation)

%%
% the system
system.sqn=[0.5 0.5];       % spin quantum numbers
system.interactions={1,0, 'z','e', nu_electron;           % Zeeman term of electron in simulation frame
                     2,0, 'z','e', nu_muon;           % Zeeman of muon
                     1,2, 'x','x', A_iso;                % isotropic HF
                     1,2, 'y','y', A_iso;                % isotropic HF
                     1,2, 'z','z', A_iso};                % isotropic HF
                     
system.init_state='ez'; % LF mode (muon in Iz eigenstate)
system.eq = 'init';  % equilibrium state is the same as initial state

% the options
options.relaxation=0;       % tells SPIDYAN whether to include relaxation (1) or not (0)
options.down_conversion=0;  % downconversion of signal (1) or not (0)
% options.det_op={'ez','ze','xe', 'ex', 'ye' }; % detection operator
options.det_op={'ez','ze',Iz_t, d34_t, d12_t , d13_t, d24_t }; % detection operator
options.det_op={d34_t, d12_t , d13_t, d24_t }; % detection operator
options.det_op={'ez',d34_t, d12_t , d13_t, d24_t }; % detection operator

options.labframe = 1;     % lab frame simulation is on
options.awg.s_rate = 500;  % can be switched on to improve accuracy

% the sequence
sequence.tp=[0, 500] ;     % vector with event lengths in ns
sequence.nu1=[0,0.95];         % amplitude, here in mT, linearly polarized
% sequence.nu1=[0,0.000001];         % amplitude, here in mT, linearly polarized
% sequence.nu1=0;         % amplitude
sequence.frq=nu_uw;  % frequency of pulse
sequence.t_rise=0;           % rise time of chirp pulses

sequence.detection=ones(1,length(sequence.tp)); % detection always on
% sequence.excite={[0 0]};      % must be commented out! otherwise it will
% not swallow the custom excitation operator

[experiment,options] = triple(sequence, options); % builds experiment (pulses)
experiment.pulse{end}.xop = exc_op;
[system, state, options]=setup(system,options); % build system 
[state, detected_signal, experiment]=homerun(state,system,experiment, options,[]);  % propagate
    
signalc=cell(1,length(B_vec)); % creates an empty cell to store results
allsignal=cell(1,length(B_vec)); % creates an empty cell to store results

tic;

act_nu_e = B0*ge;
act_nu_mu = B0*gmu;
system.interactions{1,end}= act_nu_e; 
system.interactions{2,end}= act_nu_mu;   

parfor k=1:length(B_vec)
    
    systeml=system;                          % copy system into loop
    act_nu_e = B_vec(k)*ge;
    act_nu_mu = B_vec(k)*gmu;
    systeml.interactions{1,end}= act_nu_e; 
    systeml.interactions{2,end}= act_nu_mu;   
    
    sequencel = sequence;
    
    [experimentl,optionsl] = triple(sequencel, options); % builds experiment (pulses)

    experimentl.pulse{end}.xop = exc_op;

    [systeml, state, optionsl]=setup(systeml,optionsl); % build system 
    
    [state, detected_signal, experimentl]=homerun(state,systeml,experimentl, optionsl,[]);  % propagate
    
    signalc{k}=detected_signal.sf;   % keep signals from simulation frame
    
    
    
    allsignal{k} = detected_signal;
end

toc
signal=zeros(size(signalc{1}));                  % creates empty signal
t=linspace(0,sum(experiment.tp),length(signal)); % creates time axis

% integrate
intgs = zeros(length(options.det_op),length(B0));
% filter
signalc_f=cell(1,length(B_vec)); % creates an empty cell to store results

f_cutoff_GHz = 1.0;
cutoff=f_cutoff_GHz/((1/(t(2)-t(1)))/2);

% filter is a lowpass butterworth filter of 3rd order
N=3;
[B,A] = butter(N,cutoff,'low');

for k=1:length(B0)
    signal=signal+signalc{k};      % combines all signals
    
    intgs(:,k) = sum(real(signalc{k}),2)/length(signal);
    
    % first copy over data
    signalc_f{k} = signalc{k};
    
    % then filter each signal
    for jj=1:size(signalc{k},1)
        N_pad = 100;
        tmpsig = signalc_f{k}(jj,:);
        len = length(tmpsig);
        tmpsig = [tmpsig(N_pad+1:-1:1+1) , tmpsig, tmpsig(len-N_pad+1-1:len-1)];
        tmpsig_f = filtfilt(B,A, tmpsig);
        signalc_f{k}(jj,:) = tmpsig_f(N_pad+1:len+N_pad);
    end
    
end

%%

%%
figure(1); clf; hold on
plot(t,real(signal))
% plot(t,abs(real(signal(1,:))+1i*real(signal(2,:))))
title('Signal in simulation frame')
xlabel('time [ns]')
ylabel('<I_{i}>')
legend('<I_{z}>','<S_{z}>')

offset = 0.25;
figure(2); clf; hold on
for ii=1:length(signalc)
%     plot(allsignal{ii}.t,real(signalc{ii}(4,:))+ii*offset)
    plot(allsignal{ii}.t,real(signalc{ii}(1,:))+ii*offset)
end
xlabel('time [ns]')
ylabel('<I_{i}>')
legend('<I_{z}>','<S_{z}>')


%% plot the filtered Iz operator (MuSR observable) plus the transition selective populations %Q?: filtered Iz operator?

% select dataset
trc_idx = find(B_vec == B0);

% analytical LF solution

w23 = 2*S;

LF_AC = 2*cos(alpha)^2*sin(alpha)^2*cos(2*pi*w23*t);
LF_DC = (1+(cos(alpha)^2-sin(alpha)^2)^2)/2;

figure(99); clf; hold on;
plot(t, LF_AC + LF_DC)
Iz_idx = 1;
plot(t,real(signalc{trc_idx}(Iz_idx,:)))

% plot(t,real(signalc{trc_idx}(6,:)))

% matches!



% point for data normalization : first point
normpoint = 1;


Iz_idx = 1;
sel_idx = [2, 3, 4, 5];


% Iz_norm = signalc_f{trc_idx}(Iz_idx,normpoint);
Iz_norm = LF_DC(1);
% scale the transition-wise observables by the corresponding Iz_t transition weight
d_scalings = [trace(Iz_t*d34), trace(Iz_t*d12) , trace(Iz_t*d13), trace(Iz_t*d24)];
d_norm = sum(signalc{trc_idx}(sel_idx,normpoint).*d_scalings.');

scaling = Iz_norm/d_norm; %voila, this is 1.0 by proper normalization with the trace of Iz_t

% d_scalings = d_scalings * scaling; % this scaling is no longer needed


% analytical solution for rabi oscillation frequency
nu_rabi_MHz = (cos(alpha) * abs(gmu) + sin(alpha) * abs(ge))*sequence.nu1(end)/2;

t_rabi = 1/nu_rabi_MHz*1e3;

figure(200); clf; hold on;
plot(t/t_rabi,real(signalc_f{trc_idx}(Iz_idx,:)))
ii_scaling = 1;
for ii=sel_idx
%    plot(t,real(signalc{trc_idx}(ii,:)) * d_scalings(ii_scaling))
   plot(t/t_rabi,real(signalc_f{trc_idx}(ii,:)) * d_scalings(ii_scaling))
   ii_scaling = ii_scaling + 1;
end
% plot(t,sum(real(signalc{trc_idx}(sel_idx,:).*d_scalings.'),1))
plot(t/t_rabi,sum(real(signalc_f{trc_idx}(sel_idx,:).*d_scalings.'),1))

% save('Rabi_onres','signalc_f','d_scalings', 't_rabi', 't');

%% look at ingts etc

figure(3); clf;
plot(B0*1e4, intgs)
xlabel('Field [G]')
ylabel('Time integral')
[ma, maid] = max(intgs(2,:));
legend('<Iz>','<Sz>')

p_id = maid;
% p_id = find(B_vec*1e4 > 1413 ,1);

ylims = [-1.1,1.1];
figure(4); clf;
subplot(511)
hold on
plot(t, real(real(signalc{1}(1,:))))
plot(t, real(real(signalc{p_id}(1,:))))
plot(t, real(real(signalc_f{1}(1,:))))
ylim(ylims)
ylabel('<Iz>')

subplot(512)
hold on
plot(t, real(real(signalc{1}(2,:))))
plot(t, real(real(signalc{p_id}(2,:))))
xlabel('Time [ns]')
ylabel('<Sz>')
ylim(ylims)


%% get the traces in spectral domain


tdy = zeros(length(allsignal{1}.signals.t{end}),length(signalc));
tdx = t;

for ii=1:length(signalc)
%     tdy(:,ii) = real(signalc{ii}(1,:));
%     tdy(:,ii) = allsignal{ii}.signals.dc{end}(3,:);
    % Ix
%     tdy(:,ii) = allsignal{ii}.signals.sf{end}(4,:);
%     tdy(:,ii) = tdy(:,ii) - mean(tdy(:,ii));
    
    % Iz
    tdy(:,ii) = allsignal{ii}.signals.sf{end}(1,:);
    
    
    % 34 population
    tdy(:,ii) = allsignal{ii}.signals.sf{end}(2,:);
    
end


% get spectra
[fdy, fdx] = getmultspec(tdy,tdx,2,0);

fdx = fdx;
p_idx = fdx < 100 & fdx > -100;

i2Dcut_nox(fdx(p_idx),B_vec,abs(fdy(p_idx,:)),22);
%xlim(0,100);

figure(200);
clf();
contour(fdx(p_idx),B_vec,abs(fdy(p_idx,:))',22)
% ylim([1370, 1430])
% xlim([-80, 80])
xlabel('f [GHz]')
ylabel('delay [ns]')

[~, mod_idx] = min(abs(fdx - wa*1e3));

figure(201);
clf();
plot(dly_vec, abs(fdy(mod_idx,:)))
xlabel('delay tau [ns]')
ylabel(' <Ix> at 21 MHz [a.u.]')


% get the indirect dim spectra
ind_trace = abs(fdy(mod_idx,:)).';
ind_trace = ind_trace - mean(ind_trace);
[fdy2, fdx2] = getmultspec(ind_trace,B_vec,2,0);
figure(202);
clf()
plot(fdx2, abs(fdy2))

