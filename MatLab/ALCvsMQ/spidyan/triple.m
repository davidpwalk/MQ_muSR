function [experiment,options] = triple(sequence, options)
% [experiment,options] = triple(sequence, options)
% 
% The function triple uses the structure sequence to create a pulse
% experiment that can be used with the function homerun. The 
% sweep rate of the pulse is linear, unless SPIDYAN is asked
% to compensate for the resonator profile. In that case either the
% resonator parameters or a custom resonator profile has to be given.
% If options.resonator.active is true, the resonator profile is applied to
% the calculated pulse. The pulses before convolution are stored. The
% resonator profile is also stored in the structure experiment.
% 
% A brief explanation of setting up the input can found below. For a more
% detailed instruction please refer to the documentation.
% 
% Definition of pulse sequence:
% sequence.tp       vector with event (pulses, delays, detection) lengths
% sequence.nu1      vector with pulse amplitudes
% sequence.Q        vector with critical adiabaticites for chirp pulses
% sequence.beta     vector with flip angles, in rad
% sequence.frq      cell array containing vectors with two elements,
%                   initial and final frequency of pulses
% sequence.type     cell array containing strings to identify type of
%                   requested pulse. Can be 'rectangular', 'chirp', 
%                   'linear', 'HS' or 'WURST'.
% sequence.phase    vector with phases of pulses in rad
% sequence.t_rise   rise time of the pulses in ns
% sequence.excite   cell array with boolean vectors. can be used to excite
%                   selected spins during pulses
% sequence.pcycle   cell array containing matrices for phase cycling. see
%                   below
% sequence.detection   bolean vector, elements which are to be detected
% sequence.stepwise_evolution   boolean vector, if true, stepwise
%                               propagation during free evolution events
% sequence.WURSTN   Vector with parameter for definition of WURST pulses
% sequence.HSbeta   Vectors with parameters for definition of HS pulses
% sequence.HSorder
% sequence.HSorder2
% 
% A call can be as simple as this:
% 
% sequence.tp = [20 100 20 200]
% sequence.nu1 = [40 0 40]
% sequence.frq = 1
% sequence.detection = [0 0 0 1]
% 
% This creates a two pulse experiment with two rectangular pulses at 1 GHz
% (in the simulation frame) and with a pulse amplitude of 40 MHz. The
% pulses are separated by 100 ns and can detect the signal during the last
% 200 ns period.
% The ordering in the vectors corresponds to the ordering in .tp. Pulses
% are identified through .nu1, .beta and .Q. Free evolution events must
% have 0 at their corresponding position in all three vectors. 
% 
% A more complex call could look something like this.
% 
% sequence.tp = [10 50 100 100 50 500]
% sequence.beta = pi    % equal to [pi 0 0 0 0 0]
% sequence.Q = [0 0 5]
% sequence.nu1 = [0 0 0 0 40]
% sequence.t_rise = [0 0 20 0 15]
% sequence.type = {'rectangular',[],'chirp',[],'HS'}
% sequence.frq  = {1, [], [0.8 1.2], [], [1.2 0.8]}
% sequence.HSbeta = 10       % since there is only one HS pulse we could
% sequence.HSorder = 2       % also write sequence.Hsbeta = 10
% sequence.detection = [0 0 0 0 0 1]
% sequence.stepwise_evolution(6) = 1
% sequence.pcycle{1}=[0, 1; pi, -1]
% sequence.pcycle{3}=[0, 1; pi, -1]
% 
% This will create a three pulse sequence with the first pulse being a
% rectangular pulse with a flip angle of pi, the second pulse (3rd event)
% being a chirp pulse with an adiabaticity factor of 5 and sweeping from
% 0.8 to 1.2 GHz and the third pulse being an hyperbolic secant pulse with
% a pulse amplitude of 40 MHz, a sweeprange from 1.2 to 0.8 GHz, an order
% of 2 and a beta of 5. If we had multiple pulses we could also write
% .HSbeta and HSorder in vector form, allowing to define each pulse
% individually, e.g.
% 
% sequence.HSbeta = [0 0 10 0 5]
% sequence.HSorder = [0 0 2 0 4]
% 
% The same it .frq can be simplified if all pulses have the same sweep
% widht:
% 
% sequence.frq = [0.8 1.2]
% 
% Now all pulses in the sequence are swept from 0.8 to 1.2 GHz regardless
% of their type.
% In the above example we only store expectation values during the last
% event. Also we apply phase cycling to the 1st and 3rd pulse. The first
% element per line corresponds to the phase and the second to the detection
% channel.
% By adding the line
% 
% sequence.excite = {[1 1 1] [] [1 0 0] [] [1 0 0]}
% 
% it is possible to control the excitation operator. The first vector
% means, that all three spins in the system are excited by the first pulse,
% while the other two pulses excite the first spin (as it is defined in
% system.sqn) only.
% 
% Other optional parameters can be specified with options:
% 
% Simulation is to be carried out in the labframe and the frequency for
% upconversion is 8 GHz:
% 
% options.labframe = 1
% options.LO = 8
% 
% Simulate the effect of a resonator, with center frequency 9 GHz (must 
% always be given in the lab frame) and a quality factor of 50, on the
% pulses:
% 
% options.resonator.active = 1
% options.resonator.f0 = 9
% options.resonator.Ql = 50
% 
% For simulations with a resonator, options.LO must be available.
% To account for the above resonator and compensate frequency modulated
% pulses to ensure flat excitation bands over the entire sweep width,
% additionally switch on
% 
% options.resonator.comp_bw = 1
% 
% Use complex excitation waves (real part acts on the Sx operator,
% imaginary part acts on Sy
% 
%  options.complex_excitation = 1
% 
% A useful tool for looking at your pulses is plot_pulses:
% 
% plot_pulses(experiment,options)
% 
% It will plot the polses as computated with triple. By default only the
% real part of all pulses in the time domain is depicted. Configurable
% parameters allow to look at the frequency domain instead (or time and
% frequency) and to select pulses:
% 
% options.plot_domain   'time', 'frequency' or 'both'. Selects what plots
%                       to display, default is to 'time'
% options.plot_pulse    boolean vector, ordering corresponds to
%                       experiment.tp, 1 for plotting, 0 for no plotting,
%                       if empty or not defined, all pulses are plotted.
% options.plot_imaginary_part   0/1, plot imaginary part of wave, default 
%                               is 0
% 
% 
% St. Pribitzer, 2015

% default parameters for AWG
defaultawg.s_rate=12;
defaultawg.vert_res=10;
defaultawg.channels=1;

% some defaule rise time
default_t_rise=20;

% checks the input
if ~isfield(sequence,'nu1')
    sequence.nu1=[];
end
if ~isfield(sequence,'beta')
    sequence.beta=[];
end

if isfield(options,'awg') && ~isempty(options.awg)
    awg=options.awg;
    a=0;
    % checks if all AWG fields are specified and assumes default values if
    % necessary
    if ~isfield(options.awg,'s_rate') ||  isempty(options.awg.s_rate) 
        awg.s_rate=defaultawg.s_rate;
        a=1;
    end
    if ~isfield(options.awg,'vert_res') ||  isempty(options.awg.vert_res) 
        awg.vert_res=defaultawg.vert_res;
        a=1;
    end
    
    awg.channels = 1;
      
    if a
        % stores the new AWG if fields were missing
        options.awg=awg;
        fprintf('AWG was not fully defined, assumed default values for missing fields.\n');
    end
    
else
    % If no AWG is given, assumes default values for all fields
    awg=defaultawg;
    fprintf('No AWG given, assuming default values.\n');
    options.awg=awg;
end

% sets the awg correctly for complex excitation or resonators. those
% require both AWG channels, one for the real part, the other one is the
% imaginary part
if isfield(options,'resonator') && isfield(options.resonator, 'active') && ~isempty(options.resonator.active) && (~isscalar(options.resonator.active) || options.resonator.active~=0)
    options.awg.channels=2;
end
if isfield(options,'complex_excitation') && options.complex_excitation
    options.awg.channels=2;
end

% creation of experiment
for k=1:length(sequence.tp)
    
    % checks if current event is a pulse or free evolution
    if ((isempty(sequence.nu1) || k>length(sequence.nu1) || ... 
            sequence.nu1(k)==0)) && ((isempty(sequence.beta) || ...
            k>length(sequence.beta) || sequence.beta(k)==0)) && ...
            (~isfield(sequence,'Q') || isempty(sequence.Q) || ...
            k>length(sequence.Q) || sequence.Q(k)==0)
            
        
        % for free evolution only the time is stored and the structure
        % pulse of the corresponding element is empty
        experiment.tp(k)=sequence.tp(k);
        experiment.pulse{k}=[];
        
        % creation of time step size if required
        if ~isfield(experiment,'dt') || isempty(experiment.dt)
            experiment.dt=1/options.awg.s_rate;
        end
        
        if exist('ringer','var') && ringer
            if experiment.tp(k)-offset<0
                warning('Event %1.0f: Inter-pulse delay or free evolution too short. Due to ringing the previous element had to be extened. This extension is longer than the current element. Your pulse timings might have changed. \n', k)
            else
                experiment.tp(k)=experiment.tp(k)-offset;
                ringer=0;
                fprintf('Event %1.0f: Free evolution length had to be changed from %4.2f to %4.2f ns due to ringing from the preceding pulse. \n',k,sequence.tp(k),experiment.tp(k))
            end
            
        end
    
    else
    % builds the pulse
        
        %  for a pulse SPIDYAN checks if it can find a rise time
        if isfield(sequence,'t_rise') && ~isempty(sequence.t_rise)
            if isscalar(sequence.t_rise)
                pulse.t_rise=sequence.t_rise;
            elseif k<=length(sequence.t_rise)
                pulse.t_rise=sequence.t_rise(k);
            
            else
                % if no rise time is given for the current event, SPIDYAN 
                % assumes the default rise time
                pulse.t_rise=default_t_rise;
                fprintf('Event %1.0f: No rise time given, assuming default rise time of %4.1f ns.\n', k, default_t_rise);
            end 
        else   
            % if no rise time is given, SPIDYAN assumes the default rise time
            pulse.t_rise=default_t_rise;
            fprintf('Event %1.0f: No rise time given, assuming default rise time of %4.1f ns.\n',k, default_t_rise);  
        
        end
        
                      
        if isfield(sequence,'nl_fwd') && ~isempty(sequence.nl_fwd)
            pulse.nl_fwd = sequence.nl_fwd;
        end
        
        % checks whether the correct input parameters are given for HS or
        % WURST and writes the label of the pulse type
        if isfield(sequence,'type') && ~isempty(sequence.type)
            % first checks if each pulse is already specified! If not,
            % unspecified pulses are assumed to be of the same type as the
            % first pulse
            if k <= length(sequence.type) && ~isempty(sequence.type{k})
                pulse.type = sequence.type{k};
                index = k;
            else
                pulse.type = sequence.type{1};
                fprintf('Event %1.0f: Type of pulse not specified, assuming same as first event.\n',k);
                index = 1;
            end
            
            % checks HS pulse
            if strcmp(pulse.type,'HS')
                try
                    pulse.HSbeta =sequence.HSbeta(index);
                    pulse.HSorder = sequence.HSorder(index);
                    if isfield(sequence,'HSorder2') && ~isempty(sequence.HSorder2(index))
                        pulse.HSorder2 = sequence.HSorder2(index); % only for asymmetry
                    end
                catch
                    error('Event %1.0f: Problem with HS pulse arguments.',k)
                end
                if (k<=length(sequence.nu1) && sequence.nu1(k) == 0) || (length(sequence.nu1) == 1 && sequence.nu1 == 0)
                    error('Event %1.0f: No pulse amplitude for HS pulse provided.',k)
                end
                
            % checks WURST pulse
            elseif strcmp(pulse.type,'WURST')
                try
                    pulse.WURSTN = sequence.WURSTN(index);
                catch
                    error('Event %1.0f: Problem with WURST pulse arguments. pulse.WURSTN is missing.\n',k)
                end
                if (k<=length(sequence.nu1) && sequence.nu1(k) == 0) || (length(sequence.nu1) == 1 && sequence.nu1 == 0)
                    error('Event %1.0f: No pulse amplitude for WURST pulse provided.',k)
                end
            end
            
        % if spin type is not specified at all, a simple chirp/rectangular
        % pulse is assumed. No amplitude modulation function (except for
        % smoothing flanks) is assumed.
        else
            fprintf('Event %1.0f: Type of pulse not specified, assuming linear/rectangular pulse.\n',k);
            pulse.type = 'linear';
        end
        
        % looks for pulse strength of current event
        if isfield(sequence,'nu1') && ~isempty(sequence.nu1) && ...
                k<=length(sequence.nu1) && sequence.nu1(k)~=0
            pulse.nu1=sequence.nu1(k); 
        else
            % if nu1 is not given or zero, SPIDYAN calculates the required
            % pulse strength from the flip angle
            pulse.nu1=[];
        end
        
        % looks for flip angle of current pulse
        if isfield(sequence,'beta') && ~isempty(sequence.beta) && ...
                k<=length(sequence.beta) && sequence.beta(k)~=0
            pulse.beta=sequence.beta(k);
        else 
            % if beta is not given or zero, SPIDYAN calculates the flip
            % angle from pulse length and strength
            pulse.beta=[];  
        end
        
        % In case both an amplitude and flip angle are give, a warning is
        % output. The flipangle is ignored.
        if ~isempty(pulse.nu1) && ~isempty(pulse.beta)
            fprintf('Event %1.0f: Pulse amplitude and flip angle given, using pulse amplitude.\n',k);
        end
        
        % Gets pulse length. In the case its empty or zero, the pulse
        % length will be calculated later from flip angle and amplitude
        % (only available for chirp pulses)
        if isfield(sequence,'tp') && ~isempty(sequence.tp) && ...
                k<=length(sequence.tp) && sequence.tp(k)~=0
            pulse.tp=sequence.tp(k);      
        else
            % if pulse length is not available or zero, it is calculated
            % from flip angle and pule strength
            pulse.tp=[];
        end
        
        % Next we need the excitation frequency. Three cases are supported:
        % - sequence.frq is a cell array of length 1: the cell content can
        %   be a vector or scalar. The same is applied for all pulses
        % - sequence.frq is a cell array with the same length as
        %   sequence.tp, each pulse has its own excitation band
        % - just a vector with one or two elements, same excitation
        %   frequencies are applied to all pulses.
        if isfield(sequence,'frq') && ~isempty(sequence.frq) && ...
                iscell(sequence.frq)
            
            % tests if frq is the same for every pulse
            if length(sequence.frq)==1
                
                if length(sequence.frq{1})==1 % for hard pulses 
                    pulse.frq=[sequence.frq{1} sequence.frq{1}];
                else % for chirps
                    pulse.frq=sequence.frq{1};
                end
                
            else   
            % or different for each event 
            
                if length(sequence.frq{k})==1 % for hard pulses 
                    pulse.frq=[sequence.frq{k} sequence.frq{k}];
                else % for chirps
                    pulse.frq=sequence.frq{k};
                end
            end        
            
        elseif isfield(sequence,'frq') && ~isempty(sequence.frq) && ...
                isvector(sequence.frq)
            
            % creates the excitation band in case sequence.frq is given as a
            % vector
            if length(sequence.frq)==1
                pulse.frq=[sequence.frq sequence.frq];
            else
                pulse.frq=sequence.frq;
            end
            
        else
            error('Event %1.0f: Error: No information on excitation bandwidth given.', k);
        
        end
        
        % now another check for whether a Q is available, ignored if either
        % flip angle or pulse amplitude are given.
        if isfield(sequence,'Q') && ~isempty(sequence.Q)
            if k<=length(sequence.Q)
                pulse.Q=sequence.Q(k);
            else
                pulse.Q=[];
            end
        else
            pulse.Q=[];
        end

        % looking if phase is specified for pulse
        if isfield(sequence,'phase') && ~isempty(sequence.phase)
            if length(sequence.phase)==1      
                pulse.phase=sequence.phase;
            elseif k<=length(sequence.phase)
                pulse.phase=sequence.phase(k);
            else
                pulse.phase=0;
            end            
        else
            pulse.phase=0;
        end
        
        % number of current event, required for debugging and feedback to 
        % user
        pulse.eventno=k;
        
      
        % excitation bandwidth
        dfrq=abs(pulse.frq(2)-pulse.frq(1));
              
        % in case parameters are missing, they are calculated now
        if ~isempty(pulse.tp) && isempty(pulse.nu1) &&  ~isempty(pulse.Q) && isempty(pulse.beta)
            
            if dfrq~=0 
                pulse.nu1=1000*sqrt(pulse.Q*dfrq/(2*pi*pulse.tp));
            else
                error('Event %1.0f: Can only create pulse from adiabaticity factor for a chirp.', pulse.eventno);
            end
            
        elseif ~isempty(pulse.tp) && isempty(pulse.nu1)
            
            % calculates pulse amplitude
            if ~isempty(pulse.beta)
                pulse.nu1=balk(pulse.beta,dfrq,pulse.tp, pulse.eventno);
            else
                error('Event %1.0f: Only pulse amplitude given.', pulse.eventno);
            end
        elseif isempty(pulse.tp) && ~isempty(pulse.nu1)
            
            % calculates pulse length
            if ~isempty(pulse.beta)
                pulse.tp=pickoff(pulse.beta,dfrq,pulse.nu1,pulse.eventno);
            else
                error('Event %1.0f: Only pulselength given.', pulse.eventno);
            end
            
        elseif isempty(pulse.tp) && isempty(pulse.nu1)
            
            error('Event %1.0f: Not enough pulse information given.', pulse.eventno);
        end
        
        % and now the calculation of the actual Q for chirp pulses
        if dfrq~=0
            
            % calculates adiabaticity factor for linear chirps
            pulse.Q=2*pi*(pulse.nu1/1000)^2*pulse.tp/dfrq;
            
        end
        
        pulse.nu_init=pulse.frq(1);
        pulse.nu_final=pulse.frq(2);
                
        % now that everything is setup, the waves can be computed! Repeated
        % for phasecycling
phase_offset = pulse.phase;
        if isfield(sequence,'pcycle') && k<=length(sequence.pcycle) && ~isempty(sequence.pcycle{k}) 
            
            pulse.phase=phase_offset+sequence.pcycle{k}(1,1);
                
            [wave, options]=build_pulse(pulse, options);
            
            for kk=2:size(sequence.pcycle{k},1)
                pulse.phase=phase_offset+sequence.pcycle{k}(kk,1);
                
                pulse_pcyc=build_pulse(pulse, options);
                     
                wave.binary(kk,:)= pulse_pcyc.binary;
                wave.arbitrary(kk,:)= pulse_pcyc.arbitrary;
                               
            end
            
            wave.pcycle=sequence.pcycle{k};
        else
            [wave, options]=build_pulse(pulse, options);
            wave.pcycle=0;
        end
        
        % some missing parameters are written to wave
        if dfrq==0
            wave.type = 'rectangular';
            wave.Q = [];
            wave.beta = pulse.beta;
        elseif strcmp(pulse.type,'chirp') || strcmp(pulse.type,'linear')
            wave.Q = pulse.Q;
            wave.beta = pulse.beta;
        elseif strcmp(pulse.type,'HS')
            wave.Q= [];
            wave.beta = [];
            wave.HSbeta = pulse.HSbeta;      % beta value for HS pulse
            wave.HSorder = pulse.HSorder;
            try
                wave.HSorder2=pulse.HSorder2;
            catch
                wave.HSorder2=[];
            end
        elseif strcmp(pulse.type,'WURST')
            wave.WURSTN = pulse.WURSTN;
        end
        wave.nu1=pulse.nu1;
        wave.phase = pulse.phase;
        
        % convolution of pulse with resonator profile if requested
        if isfield(options,'resonator') && isfield(options.resonator,'active') && options.resonator.active
            if exist('ringing','var') && ringer
                warning('Two pulses without delay or interpulse not sufficient to compensate for ringing. Ringing arises from convultion of the resonator profile with the original chirp. Event lengths will no longer fit and are shifted.\n')
            end
            wave=resonator(wave,options);
            ringer=1;
            offset=wave.tp-wave.tp_orig;
            if ~isfield(sequence,'tp_orig')
                experiment.tp_orig=sequence.tp;
            end
        end
        
        % and the excitation operator is written to the wave
        if isfield(sequence,'excite') && ~isempty(sequence.excite)
            
            if ~iscell(sequence.excite)
                excite=sequence.excite;
            elseif length(sequence.excite)==1  
                excite=sequence.excite{1};
            elseif k<=length(sequence.excite)
                excite=sequence.excite{k};
            else
                excite=1;
            end
            
            if size(excite,1)==size(excite,2) && size(excite,1) ~= 1
                wave.xop=excite;
            end           
        else
            excite=[];
        end
        
        wave.excite=excite;
        
        % stores the pulse into experiment
        experiment.pulse{k}=wave;
        experiment.tp(k)=experiment.pulse{k}.tp;
        experiment.dt=experiment.pulse{k}.dt;     

    end
end

% creates indices for homerun
experiment.i=1; % counter for event position
experiment.o=0; % offset of timeaxis, updated after every event in homerun

if isfield(sequence,'detection') && ~isempty(sequence.detection)
    s = zeros(1,length(experiment.tp));
    s(sequence.detection==1)= 1;
    experiment.detection = s;
else
    experiment.detection=zeros(1,length(experiment.tp));
end

if isfield(sequence,'stepwise_evolution') && ~isempty(sequence.stepwise_evolution)
    s = zeros(1,length(experiment.tp));
    s(sequence.stepwise_evolution==1)= 1;
    experiment.stepwise_evolution = s;
else
    experiment.stepwise_evolution=experiment.detection;
end

% if resonator is available it is stored to experiment as well
if isfield(options,'resonator') && isfield(options.resonator,'active') && ~isempty(options.resonator.active) && (~isscalar(options.resonator.active) || options.resonator.active~=0)
    experiment.resonator=options.resonator;
end

if isfield(sequence,'pcycle')
    experiment.pcycle=sequence.pcycle;
end

end

