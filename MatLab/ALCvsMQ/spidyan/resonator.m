function [chirp] = resonator(chirp, options)
% [chirp] = resonator(chirp, options)
% 
% The function resonator applies the resonator profile to a given pulse.
% The resonator has to be given with its quality factor and resonance
% frequency in the lab frame. If no resonator is given default values are
% assumed. The input pulse has to be complex for good results. If resonator
% is called with triple, this is already correctly set. For simulations
% carried out in the local oscillator frame, the resonator profile is
% converted. The frequency difference from the local oscillator frame to
% the lab frame can be given with options.dwnconversion, default is 8 GHz.
% The original wave before passing through the resonator is returned.
% 
% Input:
% chirp                     structure which contains all information on the
%                           waveform.
% chirp.arbitrary           waveform as double numbers, normalized to range
%                           +/-1
% chirp.binary              waveform in integer format 0...2^vert_res - 1
% 
% options       structure containing additional parameters
% resonator     structure containing the resonator settings.
%   resonator.comp_bw   if true, SPIDYAN compensates the pulse for the
%                       resonator
%   resonator.Ql        quality factor of the resonator.
%   resonator.f0        resonance frequency of the resonator.
%   resonator.s_rate     sampling rate of the resonator profile.
%   resonator.range     optional, frequency axis [GHz] for custom resonator
%                       profiles, if not given, the resonator profile
%                       is calculated from quality Ql, f0 and s_rate.
%   resonator.nu        optional, resonator profile of custom resonator.
%   .resonator.active   vector, containing indices of events. Resonator
%                       profile is applied to pulses found in .active, e.g.
%                       [1 3] will change the first and third event only.
% 
% awg           structure with information on awg, optional
%               if not given default values will be applied.
%   awg.s_rate      sampling rate [GS/s], default is 12 GS/s
%   awg.vert_res    vertical resolution [bit], default is 10 bit
%   awg.phase       phase shift in rad, optional, default is 0
% 
% LO            frequency of the local oscillator [GHz], required for
%               computations in the lab frame.
% labframe      0 or 1, default is 0. If true SPIDYAN simulates the
%               experiment in a local frame rotating at
%               options.LO+system.nu01.
% 
% Output:
% chirp                     structure which contains all information on the
%                           waveform.
% chirp.arbitrary           waveform after convolution with resonator
%                           profile
% chirp.binary              waveform in integer format 0...2^vert_res - 1
% chirp.dt                  time step [ns]
% 
% St. Pribitzer, 2014 

% optional parameter filtertruncation
if isfield(options,'filtertruncation')
    filtertrunc=option.filtertruncation;
else
    filtertrunc=1000;
end

% default resonator in case parameters are missing
defaultresonator.s_rate=36*3;
defaultresonator.Ql=60;
defaultresonator.f0=9.3;
dec_fac=0.9;

if ~isfield(options.resonator,'s_rate') || isempty(options.resonator.s_rate)
    
    options.resonator.s_rate=defaultresonator.s_rate;
    fprintf('No sampling rate for resonator given, assuming default value of %d GHz.\n', defaultresonator.s_rate);
    
end

if (~isfield(options.resonator,'Ql') || isempty(options.resonator.Ql) ...
        || ~isfield(options.resonator,'f0') || isempty(options.resonator.f0))...
        && (~isfield(options.resonator,'nu1') || ~isfield(options.resonator,'range'))
    
    Ql=defaultresonator.Ql;
    f0=defaultresonator.f0;
    s_rate=defaultresonator.s_rate;
    warning('One ore more resonator parameters are missing, assuming default values.');
    
elseif ~isfield(options.resonator,'nu1') || ~isfield(options.resonator,'range')
    
    s_rate=options.resonator.s_rate;
    Ql=options.resonator.Ql;
    f0=options.resonator.f0;
    
end

npts=2^14;

% creates frequency axis
if isfield(options.resonator,'nu1') && isfield(options.resonator,'range')
          f=options.resonator.range;
          Hid=options.resonator.nu1;
          s_rate=options.resonator.s_rate;
else
    f = transpose((-npts/2:npts/2-1)/npts*s_rate);
    
    % frequency response function of the resonator
    Hid=1./(1+1i*Ql*(f/f0-f0./f));
    
    % make it causal on our finite grid
    Hid = real(Hid) - 1i*imag(hilbert(real(Hid)));
    
end

%impulse response function
irf=ifft(ifftshift(Hid));
IRF=ifftshift(irf);


if isfield(options,'LO') && (~isfield(options,'labframe') || ~options.labframe)
    
    % calulates (downconversion) resonator in the rotating frame of the AWG
    LOf=options.LO;
        
    taxis = (0:npts-1)'/s_rate;
    t = taxis - npts/2/s_rate;

    % time domain signal of the LO
    LOt=cos(2*pi*LOf*t)-1i*sin(2*pi*LOf*t);

    IRFmixed = bsxfun(@times,LOt,IRF);
    
    awgs_rate=1/chirp.dt;

    filter_freq_n = awgs_rate/s_rate*dec_fac;
    [B,A] = butter(8,filter_freq_n,'low');
    IRFshifted = filtfilt(B,A,IRFmixed);

    IRFcut=ifftshift(IRFshifted.*chebwin(length(IRFshifted)));
    IRFcut=IRFcut(1:round(length(IRFshifted)/2)); 
       
    % hard cutting of time domain data
    ix = find(flipud(abs(IRFcut))>max(abs(IRFcut))/filtertrunc,1);
    IRFcut=IRFcut(1:end-ix+1);
    
    
    % puts resonator on same time grid as AWG
    if isfield(options,'awg') && isfield(options.awg,'s_rate') && ~isempty(options.awg.s_rate)
        awgs_rate=options.awg.s_rate;
    else
        awgs_rate=12;
    end

    IRFcorrected = IRFcut(1:s_rate/awgs_rate:end);
    
else
    
    % convolutes resonator with chirp in the labframe
    awgs_rate=1/chirp.dt;
    
    filter_freq_n = awgs_rate/s_rate*dec_fac;
    [B,A] = butter(8,filter_freq_n,'low');
    
    IRFshifted = filtfilt(B,A,IRF);
    IRFcut=ifftshift(IRF.*chebwin(length(IRFshifted)));
    IRFcut=IRFcut(1:round(length(IRFcut)/2));    
       
    ix = find(flipud(abs(IRFcut))>max(abs(IRFcut))/filtertrunc,1);
    IRFcut=IRFcut(1:end-ix+1);
    win=chebwin(length(IRFcut)*2);
    IRFcut=IRFcut.*win(length(IRFcut)+1:end);
     
    IRFcorrected = IRFcut(1:s_rate/awgs_rate:end);

end

% storing of pulse before convolution 
pulse=chirp.arbitrary; 

chirp.arbitrary_orig=chirp.arbitrary;
chirp.binary_orig=chirp.binary;
chirp.t_orig=chirp.t;
chirp.tp_orig=chirp.tp;

res=2^chirp.vert_res;

pulse_res = zeros(size(pulse,1),size(pulse,2)+length(IRFcorrected)-1);
for ii = 1:size(pulse,1)
     pulse_res(ii,:)=conv(pulse(ii,:),IRFcorrected);
end

% pulse_res=pulse_res(1:end-length(IRFcorrected)+1);

pulse_res = pulse_res / sum(abs(IRFcorrected));

real_wave=real(pulse_res);
imag_wave=imag(pulse_res);

% compute integer input values for AWG
real_binary=floor(res*(real_wave+1)/2);
real_binary(real_binary==res)=res-1;

imag_binary=floor(res*(imag_wave+1)/2);
imag_binary(imag_binary==res)=res-1;

real_wave=2*(real_binary/res-0.5);
imag_wave=2*(imag_binary/res-0.5);

% only the real part of the part is used by the wave_propagation. imaginary
% part is dropped.
chirp.binary=real_binary+1i*imag_binary;
chirp.arbitrary=real_wave+1i*imag_wave;

chirp.tp=length(pulse_res)*chirp.dt;
chirp.t=linspace(0,chirp.tp-chirp.dt,length(pulse_res));

end

