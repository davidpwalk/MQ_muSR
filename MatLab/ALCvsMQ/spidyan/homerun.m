function [state, detected_signal, experiment, options, tables, sigmas] = homerun(state, system, experiment, options, tables)
% [state, detected_signal, experiment, options, tables, sigmas] = homerun(state, system, experiment, options, tables)
%
% The function homerun is in charge of pragation of spin system from a
% given initial state through the user created pulse sequence.
% Input:
% state        density matrix of initial state, created with setup.m
% system       structure whcih contains the spin system. in specific the
%              free evolution Hamiltonian, eqilibrium state and relaxation
%              matrix gamma, usually created with setup.m
% experiment   structure which contains the entire sequence. created with
%              triple.m
% options      a structure which contains various optional parameters, see
%              below
% tables       optional, another structure which can contain propagators.
%              can be used to speed up simulations, if not available it is
%              recalculated and output
% 
% Output:
% state            density matrix, describing the state of the system after
%                  propagation
% detected_signal  structure 
%   .sf            n-dimensional vector, containing the time traces in the 
%                  simulation frame, n equals the number of detection 
%                  operators.
%   .dc            n-dimensional vector with the entire time trace of the 
%                  down converted signal, only available if 
%                  options.down_conversion is true
%   .t             vector, the time axis of the entire experiment.
%   .signals       structure which contains the individual signals of each 
%                  event
%   .signals.sf    cell, with the individual signals of each event, in 
%                  simulation frame.
%   .signals.dc    cell, containing the individual signals of each event,
%                  down converted, only available if down_conversion is 
%                  true.
%   .signals.t     cell, time axis of the individual events.
% experiment       updated version of experiment. 
%   .xop           the excitation operator in its matrix form (.xop) and
%   .xop_imaginary the complex excitation operator (.xop_imaginary) are
%                  added to each individual pulse (experimen.pulse{k})
%                  
% options          updated structure options
%   .down_conversion            is set to zero if not defined
%   .store_density_matrix       is set to zero if not defined
% 
% tables        structure containing lookup tables
% .U            table with precalculated Propagators in the Hilbert space
% .sigss        table with steadystate density operators.can only be used, 
%               for pulses with the same nu1 and system Hamiltonian.
% .L            table with Liouvillians.can only be used, if the exact same
%               pulse and system are used again.
% .nu1          pulse amplitude of the pulse for whiche tables are given. 
%               If not the same as the applied pulse, tables will be 
%               calculated new
% .ham          Hamiltonian of the system, used for checking if the spin
%               system has changed.
% 
% Configureable parameters:
% experiment.detection          boolean vector, position corresponds to
%                               elements in experiment.tp, 1 for
%                               computation of expectation values, 0 off
% experiment.stepwise_evolution boolean vector, position corresponds to
%                               elements in experiment.tp, only available
%                               for free evolution events, if true
%                               stepwise propagation on the grid of the AWG
%                               else propagation of state in one large time
%                               step (faster).
% 
% options.relaxation            1/0
% options.complex_excitation    1/0
% options.store_density_matrix  1/0, if true density matrices are stored
%                               for each time step, requires
%                               experiment.detection to be true for current
%                               event, can be plotted with barcscroll
% options.down_conversion        1/0 if true, the individual and the entire
%                               time trace will also be returned after down
%                               conversion. Requires a downconversion
%                               frequency system.dc_freq. If non available
%                               SPIDYAN tries to guess frequency for down
%                               conversion by taking the highest frequency
%                               in the Hamiltonian table (see setup.m)
% 
% system.dc_freq                optional, frequency that can be used for
%                               down conversion. If not available the
%                               highest frequency in the Hamiltonian tbale
%                               is used instead.
% 
% St. Pribitzer, 2015


% switch down conversion to false, if not defined
if ~isfield(options,'down_conversion') || isempty(options.down_conversion)
    options.down_conversion=0;
end

% if specific freequency for down conversion is given, takes it from
% system, else the highest frequency (assumed to be a Zeeman term) is taken
if isfield(system,'dc_freq') && ~isempty(system.dc_freq)
    dc_freq=system.dc_freq;
else
    dc_freq=system.highest_freq;
end


% sets up cells for the output
if isempty(options.detect)
    detected_signal=[];
end
signal=cell(1,length(experiment.tp));
signalf=cell(1,length(experiment.tp));
times=cell(1,length(experiment.tp));
tablesnew=cell(1,length(experiment.tp));

if isfield(options,'store_density_matrix') && ~isempty(options.store_density_matrix) && options.store_density_matrix
    sigmas = cell(1,length(experiment.tp));
else
    sigmas = [];
    options.store_density_matrix  = 0;
end

% this creates the excitation operators for the pulses, by calling a local
% function and stores it in experiment.pulse{k}.xop
[experiment, system, options] = mk_excitation_operators(experiment, system, options);

% due to the choice of the time grid, the initial state has to be rotated
% back by half a timestep
backrot = expm(1i*system.ham*experiment.dt/2);
state = backrot*state*backrot';

if nargin<4
    tables = [];
end

% create some indices used w=by wave_propagation
c=0;
r=1;

% calls the function wave propagation to propagate the system, with table
% if possible
for k=1:length(experiment.tp)
    
    if ~isempty(tables) && k<=length(tables)
        table = tables{k};
    else
        table = [];
    end 
    
    % porpagate wave
    [state, signal{k}, experiment, table, sigmas{k}] = wave_propagation(state, system, experiment, options, table);
    
    % measure length of returned
    x = size(signal{k},2);
    
    % down conversion makes sense only if the signal exceeds a certain
    % length
    if options.down_conversion && x>10
        [signalf{k},options.filter]=strike(signal{k}, experiment.t{k}, dc_freq, options);
    else
        signalf{k}=[];
    end
    
    % updating the timeaxis
    c=c+x;
    times{k}=experiment.t{k};
    r=r+1;
    
    % storing new table
    tablesnew{k}=table;
    if ~isempty(tablesnew{k})
        tablesnew{k}.ham=system.ham;
    end
       
end

% if the new tables differ from the ones provided (if any), they are stored
% and returned
if ~isequal(tablesnew,table)
        tables=tablesnew;
end

% this rearranges the detected signals and time axis and combines them, so
% that it is possible to plot them all as one
if ~isempty(options.detect) 
    % setting correct sizer of vector
    detected_signal.sf=zeros(length(options.detect),c-experiment.i+2);
    experiment.taxis=zeros(1,c-experiment.i+2);
    
    %indices again
    a=2;
    first = 1;
    
    for k=1:length(signal)
        
        % check if first signal is not empty
        if first && ~isempty(signal{k})
            
            % store first point of the first signal
            detected_signal.sf(:,1)=signal{k}(:,1);
            experiment.taxis(1,1)=times{k}(1);
            
            % now add all the others timepoints
            [~,y]= size(signal{k});
            b=a+y-2;
            detected_signal.sf(:,a:b)=signal{k}(:,2:end);
            experiment.taxis(1,a:b)=times{k}(2:end);
            a=b+1;
            first = 0;
            
        elseif ~isempty(signal{k})
            % adding other signals, this confusing index is necessary in 
            % order to avoid double counting last point of a signal and the
            % first point of the succiding signal
            [~,y]= size(signal{k});
            b=a+y-2;
            detected_signal.sf(:,a:b)=signal{k}(:,2:end);
            experiment.taxis(1,a:b)=times{k}(2:end);
            a=b+1;
        end
               
    end
    
    % down conversion of the entire signal
    [~, x]= size(detected_signal.sf);     % check if signal length is sufficient for down_conversion           
    
    if options.down_conversion && x>10
        [detected_signal.dc,options.filter]=strike(detected_signal.sf, experiment.taxis, dc_freq, options);
        signals.dc=signalf;
    else
        detected_signal.dc=[];
        signals.dc=[];
    end
    
    % and storing everything
    signals.sf=signal;
    detected_signal.signals=signals;
    detected_signal.t=experiment.taxis;
    detected_signal.signals.t=times;
       
end

end


function [experiment, system, options] = mk_excitation_operators(experiment, system, options)

% building the excitation operator
for k = 1: length(experiment.pulse)
    if ~isempty(experiment.pulse{k})
        if ((~isfield(experiment.pulse{k},'excite') || isempty(experiment.pulse{k}.excite)) && (~isfield(experiment.pulse{k},'xop') || isempty(experiment.pulse{k}.xop))) || length(experiment.pulse{k}.excite)==1
            xop = 0;
            for l = 1: length(system.sqn)
                op=1;
                for ll = 1:length(system.sqn)
                    if ll~=l
                        op = kron(op,spops(system.sqn(ll),'e'));
                    elseif ll==l
                        op = kron(op,spops(system.sqn(ll),'x'));
                    end
                end
            xop =  xop+op;    
            end
            experiment.pulse{k}.xop = xop;
        elseif isfield(experiment.pulse{k},'excite') && ~isempty(experiment.pulse{k}.excite) && length(experiment.pulse{k}.excite)~=1 && size(experiment.pulse{k}.excite,1)==1
            xop = 0;
            for l = 1: length(system.sqn)
                if l <= length(experiment.pulse{k}.excite) && experiment.pulse{k}.excite(l) == 1
                    op=1;
                    for ll = 1:length(system.sqn)
                        if ll == l
                            op = kron(op,spops(system.sqn(ll),'x'));
                        else
                            op = kron(op,spops(system.sqn(ll),'e'));
                        end
                    end
                    xop =  xop+op;
                end
            end
            experiment.pulse{k}.xop = xop;
        elseif isfield(experiment.pulse{k},'excite') && ~isempty(experiment.pulse{k}.excite) && size(experiment.pulse{k}.excite,1)==size(experiment.pulse{k}.excite,2)
        else
            if (~isfield(experiment.pulse{k},'xop') || isempty(experiment.pulse{k}.xop))
                error('Whoops, something went wrong with the excitation operator, please check.')
            end
        end
    end
end

% and here the complex_excitation operator, if requested
if isfield(options,'complex_excitation') && options.complex_excitation
    for k = 1: length(experiment.pulse)
        if ~isempty(experiment.pulse{k})
            if ((~isfield(experiment.pulse{k},'excite') || isempty(experiment.pulse{k}.excite)) && (~isfield(experiment.pulse{k},'xop') || isempty(experiment.pulse{k}.xop))) || length(experiment.pulse{k}.excite)==1
                xop = 0;
                for l = 1: length(system.sqn)
                    op=1;
                    for ll = 1:length(system.sqn)
                        if ll~=l
                            op = kron(op,spops(system.sqn(ll),'e'));
                        elseif ll==l
                            op = kron(op,spops(system.sqn(ll),'y'));
                        end
                    end
                    xop =  xop+op;
                end
                experiment.pulse{k}.xop_complex = xop;
            elseif isfield(experiment.pulse{k},'excite') && ~isempty(experiment.pulse{k}.excite) && length(experiment.pulse{k}.excite)~=1
                xop = 0;
                for l = 1: length(system.sqn)
                    if l <= length(experiment.pulse{k}.excite) && experiment.pulse{k}.excite(l) == 1
                        op=1;
                        for ll = 1:length(system.sqn)
                            if ll == l
                                op = kron(op,spops(system.sqn(ll),'y'));
                            else
                                op = kron(op,spops(system.sqn(ll),'e'));
                            end
                        end
                        xop =  xop+op;
                    end
                end
                experiment.pulse{k}.xop_complex = xop;
            else
                error('Whoops, something went wrong with the complex excitation operator, please check.')
            end
        end
    end
end

end


