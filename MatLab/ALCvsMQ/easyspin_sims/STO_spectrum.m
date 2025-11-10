% spectra for muonium in STO

clear all

% Zeeman
ge = 28.02495;

% magnetic field
B0 = 0.27;

% microwave frequency
freq = ge*B0;
freq = 3.456;
% freq = 345.6;

% orientation
crystal = [0.0, 0.0, 0.0];


% value for free electron, scaled down by three to make 1H equal to muon:
STO.g = [2.0 / 3.0];
    
STO.Nucs = '1H';

A_iso = 1.4;
A_dip = 15.4;
STO.A = [-A_dip+A_iso -A_dip+A_iso 2*A_dip+A_iso];   % MHz
% Zaher hamiltonion
% STO.A = [1.4 6.7 11.5];   % MHz


% STO.lw = 2.5; % 

Exp.mwFreq = freq;    % GHz
Exp.Harmonic = 0;    % no field modulation!
% Exp.Range = [360 380];   % mT

% Exp.nPoints = 561; % equivalent to 5 MHz frequency stepping
% Exp.nPoints = 701; % equivalent to 4 MHz frequency stepping
% Exp.MolFrame = crystal;

Opt.separate = 'transitions';

[B, spc, info] = pepper(STO, Exp, Opt);

% spc = spc / trapz(B,spc);



figure(1); clf; hold on;
plot(B, spc)


xlabel('Field [mT]')
ylabel('Wgt')

% ylim([-0.03,0.5])
% opts.vstretch = 0.65;
% opts.hstretch = 1.0;
% publishish(gca,gcf,opts)
% update_fonts2(gca,6,'Times', 6)



