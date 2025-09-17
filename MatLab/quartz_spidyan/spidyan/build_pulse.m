function [wave, options]=build_pulse(pulse,options)
% [wave, options]=build_pulse(pulse,options)
%
% Generates a pulse waveform for spin dynamics simulation with SPIDYAN. 
% Currently the following pulses are available:
%   - rectangular pulse
%   - linear chirp
%   - hyperbolic secant (HS) pulse
%   - WURST pulse
% 
% The function is usually called by triple.m but can also be used by its
% own. The input requires the two structures pulse and options.
% 
% pulse.tp          pulse length in ns
% pulse.nu_init     initial frequency of the pulse in the simulation frame
%                   in  GHz.
% pulse.nu_final    final frequency of the pulse, in GHz, can be the same
%                   as pulse.nu_init for a rectangular pulse.
% pulse.type        string, can be 'rectangular', 'chirp', 'linear' (the
%                   same as chirp), 'HS' or 'WURST'
% 
% options.awg       structure with the AWG parameters, such as sampling
%                   rate (.s_rate), vertical resolution (.vert_res) and
%                   channels (.channels)
% 
% optional Input:
% Pulse parameters for HS pulses:    pulse.HSbeta
%                                    pulse.HSorder
% Pulse parameters for WURST pulses: pulse.WURSTN
% 
% Simulation in labframe:  options.labframe
%                          options.LO
% 
% Resonator:               options.resonator.active
%                          options.resonator.f0
%                          options.resonator.Ql
%                          options.resonator.bw_compensation
% 
% For mor details on how to set up input, please refer the to manual or the
% documentation of triple.m.
% 
% Stephan Pribitzer, 2015

% checks 
if ~isfield(pulse,'scale'),
    pulse.scale=1;
end

if ~isfield(pulse,'HSorder2') && isfield(pulse,'HSorder'),
    pulse.HSorder2=pulse.HSorder;
end

awg=options.awg;

% converts vertical resolution in steps
res=2^awg.vert_res;

% converts the excitation band to the labframe if requested
if isfield(options,'labframe') && ~isempty(options.labframe) && options.labframe && isfield(options,'LO') && ~isempty(options.LO)
    pulse.nu_final=pulse.nu_final+options.LO;
    pulse.nu_init=pulse.nu_init+options.LO;
end

% creates higher sampling rate if given sampling rates is not sufficient
if awg.s_rate<2*pulse.nu_final || awg.s_rate<2*pulse.nu_init
    if pulse.nu_final>pulse.nu_init
        dt=1/(2.5*(pulse.nu_final));
    else
        dt=1/(2.5*(pulse.nu_init));
    end
    warning('Event %1.0f: Sampling rate of AWG not sufficient for requested pulse. Assuming sampling rate of %3.1f GS/s.',pulse.eventno,1/dt);
    options.awg.s_rate=1/dt;
    options.awg.changed=1;
else
    dt=1/awg.s_rate;
end

% comparison of pulse length and rise time, and adaption of latter if
% necessary
if pulse.t_rise==pulse.tp
    pulse.t_rise=pulse.tp-pulse.tp/10;
    warning('Event %1.0f: Rise time of pulse is the same as pulse length. New rise time is %3.1f ns.',pulse.eventno,pulse.t_rise);
elseif pulse.t_rise>pulse.tp
    pulse.t_rise=pulse.tp-pulse.tp/10;
    warning('Event %1.0f: Rise time of pulse is larger than pulse length. New rise time is %3.1f ns.',pulse.eventno,pulse.t_rise);
elseif pulse.t_rise*2>=pulse.tp
    warning('Event %1.0f: Rise time is half or more than the pulse length, you might want to consider changing it.',pulse.eventno);
end

% creates time axis
n=ceil(pulse.tp/dt);
wave.tp=n*dt;
wave.dt=dt;
t=0:wave.dt:wave.tp-wave.dt;
wave.t=t;

% if resonator parameters or a resonator profile are available and
% compensation is turned on, the pulses can be adapted in such a way, that
% a flat excitation profile is provided. Caution: The maximal amplitude and
% (for chirps therefore the Q) no longers necessarily correspond to the initially
% provided values. A feature for that will be added in the future.
if isfield(options,'resonator') && isfield(options.resonator,'comp_bw') && ~isempty(options.resonator.comp_bw) && options.resonator.comp_bw && isfield(options.resonator,'comp_bw') && options.resonator.active
    %check for completness
    if ~isfield(options.resonator,'s_rate') || isempty(options.resonator.s_rate)
        options.resonator.s_rate=36*3;
    end
    
    % checks if all necessary fields for the resonator are available and
    % uses measured profile, or computes profile if necessary
    if isfield(options.resonator,'nu1') && isfield(options.resonator,'range')
        resntor.range=options.resonator.range;
        resntor.nu1=options.resonator.nu1;  
    elseif isfield(options.resonator,'Ql') && isfield(options.resonator,'f0')
        npts=2^14;
        f = transpose((-npts/2:npts/2-1)/npts*options.resonator.s_rate);
        Hid=1./(1+1i*options.resonator.Ql*(f/options.resonator.f0-options.resonator.f0./f));
        if ~isfield(options,'labframe') || ~options.labframe
            resntor.range=f-options.LO;
        else
            resntor.range=f;
        end
        resntor.nu1=Hid;
        wave.resonator = resntor;
    else
        warning('Resonator compensation requested, but resonator not defined, no compensation.');
        resntor.range = [0 20];
        resntor.nu1 = [1 1];
    end
    
    % check the resonator range
    if pulse.nu_init<min(resntor.range),
        warning('Requested min. frq. of %4.2f GHz is less than min. frq. of %4.2f GHz of the known resonator profile.',pulse.nu_init,min(resntor.range));
    end
    if pulse.nu_final>max(resntor.range),
        warning('Requested max. frq. of %4.2f GHz is more than max. frq. of %4.2f GHz of the known resonator profile.',pulse.nu_final,max(resntor.range));
    end
    
    % resonator compensation for chirps and WURST
    if ~isfield(pulse,'type') || isempty(pulse.type) || strcmp(pulse.type,'chirp') || strcmp(pulse.type,'WURST') || strcmp(pulse.type,'linear')
        if ~isfield(pulse,'type') && strcmp(pulse.type,'WURST')
            warning('Event %1.0f: Bandwidth compensation for WURST type pulses not available, assuming chirp pulse instead.',pulse.eventno)   
        end
        
        pulse.type = 'chirp';
        frange=linspace(pulse.nu_init,pulse.nu_final,length(t));   % new
%         frange=pulse.nu_init:sign(pulse.nu_final-pulse.nu_init)*0.001:pulse.nu_final;
        if length(frange) < 2
            frange = ones(2,1)*pulse.nu_init;
        end
        sprofile=interp1(resntor.range,abs(resntor.nu1),frange,'pchip');
        mean_alpha=trapz(1./sprofile.^2)/t(end); % 2pi/Qref
        dt_df = 1./sprofile.^2/mean_alpha;
        wave.qref = 2*pi/mean_alpha/0.001; % division by 0.001 (trapz step size)
        
        % check for resonator nu1 in MHz
        if max(abs(resntor.nu1)) > 1
            wave.qref = wave.qref * 0.001^2; % make MHz to GHz
        end
        % get time tau vs freq from dt/df
        tau = cumtrapz(dt_df);
        
        % reparametrize tau(nu) to nu(tau)
        nu = interp1(tau,frange,t,'pchip');
        wave.nu = nu;
        
        % smoothing of the edges with a quarter cosine, the amplitude
        % modulation function
        smoothing = ones(1,length(nu));
        nr=round(pulse.t_rise/wave.dt);
        if nr>0,
            dr=nr*wave.dt;
            tr=dr:-wave.dt:0;
            edge=cos(pi*tr/(2*dr));
            smoothing(1:nr+1)=edge.*smoothing(1:nr+1);
            edge=fliplr(edge);
            smoothing(end-nr:end)=edge.*smoothing(end-nr:end);
        end

        wave.AM = interp1(frange,smoothing,wave.nu,'pchip');
        wave.FM = wave.nu;
        % calculate the wave
        real_wave=pulse.scale.*cos(2*pi*cumtrapz(nu)*wave.dt+pulse.phase);
        imag_wave=pulse.scale.*sin(2*pi*cumtrapz(nu)*wave.dt+pulse.phase);
        
    % resonator compensation for HS pulses
    elseif strcmp(pulse.type,'HS')
        
        deltaf = pulse.nu_final - pulse.nu_init;
        tcent = pulse.tp / 2;
        
        % incorporate AM
        beta = pulse.HSbeta;
        order = pulse.HSorder;
        beta_exp = log(beta*(0.5)^(1-order))/log(beta);
        cut = round(length(t)/2);
        AM(1:cut) = sech(beta^beta_exp *((t(1:cut)-tcent)/pulse.tp).^order); % sech
        order = pulse.HSorder2;
        beta_exp = log(beta*(0.5)^(1-order))/log(beta);
        AM(cut+1:length(t)) = sech(beta^beta_exp * ((t(cut+1:end)-tcent)/pulse.tp).^order);
        FM = deltaf*cumtrapz(t,AM.^2)/trapz(t,AM.^2) + pulse.nu_init;
        if isfield(pulse,'nl_fwd')
            AM = polyval(pulse.nl_fwd,AM);
        end
        frange = FM;
        smoothing = AM;
        
        sprofile=interp1(resntor.range,abs(resntor.nu1),frange,'pchip');
%         sprofile_orig = sprofile;
        
        % multiplicative
        sprofile = sprofile .* smoothing;
        
        mean_alpha=trapz(frange,1./sprofile.^2)/t(end); % 2pi/Qref
        dt_df = 1./sprofile.^2/mean_alpha;
%         wave.qref = 2*pi/mean_alpha/0.001; % division by 0.001 (trapz step size)
        wave.qref = 2*pi/mean_alpha; % division by 0.001 (trapz step size)
        
        % check for resonator nu1 in MHz
        if max(abs(resntor.nu1)) > 1
            wave.qref = wave.qref * 0.001^2; % make MHz to GHz
        end
        % get time tau vs freq from dt/df
        tau = cumtrapz(frange,dt_df);
        
        % reparametrize tau(nu) to nu(tau)
        nu = interp1(tau,frange,t,'pchip');
        wave.nu = nu;
        
        % bring back the AM to time domain
        wave.AM = interp1(frange,smoothing,wave.nu,'pchip');
        wave.FM = wave.nu;
        
        
        % calculate the wave
        real_wave=pulse.scale.*wave.AM.*cos(2*pi*cumtrapz(wave.FM)*wave.dt+pulse.phase);
        imag_wave=pulse.scale.*wave.AM.*sin(2*pi*cumtrapz(wave.FM)*wave.dt+pulse.phase);
        
    end
    
    
% computation of "pure" wave forms
elseif isfield(pulse,'type') && strcmp(pulse.type,'HS')
    % implement HS pulse with parameters
    
    deltaf = pulse.nu_final - pulse.nu_init;
  
    tcent = pulse.tp / 2; % centertime
     
    beta = pulse.HSbeta;
    beta_exp = log(beta*(0.5)^(1-pulse.HSorder))/log(beta);
    wave.AM = sech(beta^beta_exp * ((t-tcent)/pulse.tp).^pulse.HSorder); % Amplitude modulation of the pulse
    wave.FM = deltaf*cumtrapz(t,wave.AM.^2)/trapz(t,wave.AM.^2) + pulse.nu_init;  % Frequency modulation of pulse
    
    if isfield(pulse,'nl_fwd')
        wave.AM = polyval(pulse.nl_fwd,wave.AM);
    end
    
    % calculate the wave
    real_wave=pulse.scale.*wave.AM.*cos(2*pi*cumtrapz(wave.FM)*wave.dt+pulse.phase);
    imag_wave=pulse.scale.*wave.AM.*sin(2*pi*cumtrapz(wave.FM)*wave.dt+pulse.phase);
    
elseif isfield(pulse,'type') && strcmp(pulse.type,'WURST')
    % implement WURST pulse with parameters
    
    deltaf = pulse.nu_final - pulse.nu_init;
    
    tcent = pulse.tp / 2; % centertime
    
    wave.AM = 1-abs(sin((pi*(t-tcent))/pulse.tp)).^pulse.WURSTN; % WURST
    wave.FM = pulse.nu_init + deltaf/pulse.tp * t;
    
    if isfield(pulse,'nl_fwd')
        wave.AM = polyval(pulse.nl_fwd,wave.AM);
    end
    
    % calculate the wave
    real_wave=pulse.scale.*wave.AM.*cos(2*pi*cumtrapz(wave.FM)*wave.dt+pulse.phase);
    imag_wave=pulse.scale.*wave.AM.*sin(2*pi*cumtrapz(wave.FM)*wave.dt+pulse.phase);
    
else
    % implement linear chirp or rectangular pulse
    
    alpha=(pulse.nu_final-pulse.nu_init)/pulse.tp;
    
    wave.FM = pulse.nu_init+alpha*t;
    
    % smoothing of corners, corresponds to amplitude modulation
    smoothing = ones(1,length(t));
    nr=round(pulse.t_rise/wave.dt);
    if nr>0,
        dr=nr*wave.dt;
        tr=dr:-wave.dt:0;
        edge=cos(pi*tr/(2*dr));
        smoothing(1:nr+1)=edge.*smoothing(1:nr+1);
        edge=fliplr(edge);
        smoothing(end-nr:end)=edge.*smoothing(end-nr:end);
    end
    wave.AM = smoothing;
    
    % compute pulse
    real_wave=wave.AM.*cos(2*pi*(pulse.nu_init*t+(alpha/2)*t.^2)+pulse.phase);
    imag_wave=wave.AM.*sin(2*pi*(pulse.nu_init*t+(alpha/2)*t.^2)+pulse.phase);
    
end

% Can store original wave form, before discretization, for checking whether high frequency artifats are a problem
% wave.o = real_wave;

% compute integer input values for AWG
real_binary=floor(res*(real_wave+1)/2);
real_binary(real_binary==res)=res-1;

imag_binary=floor(res*(imag_wave+1)/2);
imag_binary(imag_binary==res)=res-1;

% stores complex part of the wave if requested
if awg.channels==2,
    wave.binary=real_binary+1i*imag_binary;
else
    wave.binary=real_binary;
end

% compute normalized with with vertical resolution of AWG
real_wave=2*(real_binary/res-0.5);
imag_wave=2*(imag_binary/res-0.5);

if awg.channels==2,
    wave.arbitrary=real_wave+1i*imag_wave;
else
    wave.arbitrary=real_wave;
end

% storing of some identifier fields
wave.vert_res = awg.vert_res;
wave.nu_init=pulse.nu_init;
wave.nu_final=pulse.nu_final;
wave.type = pulse.type;

end