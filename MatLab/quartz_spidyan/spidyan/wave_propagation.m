function [state, detected_signal, experiment, tables, sigmas]=wave_propagation(state, system, experiment, options, tables)
% [state, detected_signal, experiment, tables, sigmas]=wave_propagation(state, system, experiment, options, tables)
%
% The function homerun propagates the spin density matrix of a given spin
% system for a pulse or a free evolution event. If detection operators are
% defined, see setup.m, expectation values can be computed for every time
% step, detected_signal. It is also possible to store density matrices for
% each time step.
% Homerun also supports phasecycling, as defined with triple.m and complex
% excitation wave forms.
% Also output are an updated version of the structure experiment and
% tables, which contains propagators (for no relaxation) or Liovillians and
% steady state density matrices (for relaxation) and can be used for
% consecutive elements.
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
%
% INFO: It is not necessary to use the built-in MATLAB function to compute
% the trace for expecation values (trace(Sz*sigma)). This program uses the
% following method instead
%                           sum(sum(Sz.*sigma.'))
% Compared to the built-in function, this provides a big performance boost.
% Thanks go to Nino Wili for pointing this out.
%
% St. Pribitzer, 2015


if ~isempty(experiment.pulse{experiment.i}) % this checks if the current element in the sequence is a pulse
    
    wave=experiment.pulse{experiment.i};
    waves.arbitrary=wave.arbitrary;
    waves.binary=wave.binary;
    
    state0=state;
    
    %   creates structure for phase cycling
    if size(wave.pcycle,1)>1
        loopstate=zeros(size(state));
        norm=sumabs(wave.pcycle(:,2));
    end
    
    %  now this runs over the length of the phase cycle (1 if no phacecycling)
    for rr=1:size(wave.pcycle,1)
        
        state=state0;
        
        %  extracts imaginary part of the wave if .complex_ecitation is true
        if isfield(options,'complex_excitation') && options.complex_excitation
            wave.arbitrary_imag=imag(waves.arbitrary(rr,:));
            wave.binary_imag=imag(waves.binary(rr,:));
            complex=1;
        else
            complex=0;
        end
        
        wave.arbitrary=real(waves.arbitrary(rr,:));
        wave.binary=real(waves.binary(rr,:));
        
        MHz2GHz=1e-3;
        
        % create time axis
        npoints=ceil(wave.tp/wave.dt);
        t=linspace(0,npoints*wave.dt,npoints+1);
        
        if ~isfield(wave,'nu1'),
            error('Error: wave_propagation called without nu1 for the pulse.');
        end
        
        % proper amplitude scaling
        wave.arbitrary=wave.arbitrary*2*2*pi*wave.nu1*MHz2GHz; % factor 2 required because of linearly polarized irradiation
        
        % moves the wave in binary format to -vert_res/2...vert_res/2-1 for proper
        % scaling
        vert_res=2^wave.vert_res;
        binary=wave.binary-vert_res/2;
        scale=max(wave.arbitrary)/max(binary);
        
        % set up denisty matrices and free evolution Hamiltion
        sig=state;
        ham0=system.ham;
        
        if options.relaxation
            sig_eq=reshape(system.eq,size(system.eq,1)^2,1);
        end
        
        
        % get excitation operator (created with triple.m)
        xop=wave.xop;
        if isfield(options,'complex_excitation') && options.complex_excitation
            xop_complex=wave.xop_complex;
        end
        
        % loads propagators and other precomputed quantities, if available,
        % checks if pulse amplitude and system are consistent
        % else tables are computed
        % skipped for complex wave form
        if exist('tables','var') && ~isempty(tables) && tables.nu1==wave.nu1 && isequal(tables.ham,system.ham) && ~complex
            
            if options.relaxation == 0
                Utable=tables.U;
            else
                Ltable=tables.L;
                sigsstable=tables.sigss;
            end
            
        elseif options.relaxation == 0 && ~complex
            % if phasecycling is active, all possible propagators are
            % precomputed
            if size(wave.pcycle,1)>1
                for k=0:vert_res-1
                    ham1=scale*(k-vert_res/2)*xop;
                    ham=ham0+ham1;
                    U=expm(-1i*ham*wave.dt);
                    Utable{k+1}=U;
                end
                
            else % for no phase cycling it is sufficient to compute the actual number of propagators required
                for k=min(wave.binary):max(wave.binary)
                    ham1=scale*(k-vert_res/2)*xop;
                    ham=ham0+ham1;
                    U=expm(-1i*ham*wave.dt);
                    Utable{k+1}=U;
                end
            end
            
            % stores the table
            tables.nu1=wave.nu1;
            tables.U=Utable;
            tables.ham=system.ham;
            
        elseif options.relaxation == 1 && ~complex
            
            [n,~]=size(sig);
            sigsstable=cell(1,vert_res);
            Ltable=cell(1,vert_res);
            
            % if phasecycling is active, all possible propagators are
            % precomputed
            if size(wave.pcycle,1)>1
                for k=0:vert_res-1
                    ham1=scale*(k-vert_res/2)*xop;
                    ham=ham0+ham1;
                    hamm=kron(eye(n,n),ham)-kron(ham.',eye(n,n));
                    L=-1i*hamm-system.gamma;
                    sigss=system.gamma*sig_eq; % steady state solutions for the denisty matrices
                    sigss=-L\sigss;
                    sigsstable{k+1}=sigss;
                    L=expm(L*wave.dt); %calculations of the Liouvillians
                    Ltable{k+1}=L;
                end
            else % for no phase cycling it is sufficient to compute the actual number of propagators required
                for k=min(wave.binary):max(wave.binary)
                    ham1=scale*(k-vert_res/2)*xop;
                    ham=ham0+ham1;
                    hamm=kron(eye(n,n),ham)-kron(ham.',eye(n,n));
                    L=-1i*hamm-system.gamma;
                    sigss=system.gamma*sig_eq; % steady state solutions for the denisty matrices
                    sigss=-L\sigss;
                    sigsstable{k+1}=sigss;
                    L=expm(L*wave.dt); %calculations of the Liouvillians
                    Ltable{k+1}=L;
                end
            end
            
            % stores the table
            tables.nu1=wave.nu1;
            tables.L=Ltable;
            tables.ham=system.ham;
            tables.sigss=sigsstable;
            
        end
        
        % now propagation starts
        % first we need to figure out, whether computation of expextation
        % values is requested for the current element in the pulse sequence
        if ~isempty(options.detect) && experiment.detection(experiment.i)
            % if detection is requested, the output is set up
            
            ld=length(options.detect);
            detected_signal=zeros(ld,length(wave.arbitrary));
            norms=zeros(1,length(options.detect));
            
            if isfield(options,'store_density_matrix') && ~isempty(options.store_density_matrix) && options.store_density_matrix
                sigmas=cell(1,length(wave.arbitrary)+1);
                sigmas{1}=sig;
            end
            
            % computation of norms for detection
            for kk=1:length(options.detect),
                norms(kk) = sum(sum(options.detect{kk}.*options.detect{kk}));
                detected_signal(kk,1) = sum(sum(options.detect{kk}.*sig.'))/norms(kk);
            end
        end
        
        if options.relaxation==0
            % propagates without relaxation in Hilbert space
            for k=1:length(wave.arbitrary)
                if ~complex
                    % gets propagators from table
                    U=Utable{wave.binary(k)+1};
                    sig=U*sig*U';
                else
                    % if complex excitation is requested, usage of tables
                    % is not feasible, and propagators are computed for
                    % each time step
                    ham1=scale/2*(wave.binary(k)-vert_res/2)*xop+scale/2*(wave.binary_imag(k)-vert_res/2)*xop_complex;
                    ham=ham0+ham1;
                    U=expm(-1i*ham*wave.dt);
                    sig=U*sig*U';
                end
                
                if ~isempty(options.detect) && experiment.detection(experiment.i)
                    % computes expectation values
                    for j=1:ld
                        detected_signal(j,k+1)=sum(sum(options.detect{j}.*sig.'))/norms(j);
                    end
                    
                    % and stores density matrices if requested
                    if isfield(options,'store_density_matrix') && ~isempty(options.store_density_matrix) && options.store_density_matrix
                        sigmas{k+1}=sig;
                    end
                end
            end
            
        else
            % propagation with relaxation in Liouville space
            n = size(sig,1);
            sigvec=reshape(sig,n*n,1);
            
            for k=1:length(wave.arbitrary)
                if ~complex
                    % gets Liouvillians from table
                    L=Ltable{wave.binary(k)+1};
                    sigss=sigsstable{wave.binary(k)+1};
                    sigvec=sigss+L*(sigvec-sigss);
                else
                    % if complex excitation is requested, usage of tables
                    % is not feasible, and Liouvillians and state state
                    % density matrices are computed for each time step
                    ham1=scale/2*(wave.binary(k)-vert_res/2)*xop+scale/2*(wave.binary_imag(k)-vert_res/2)*xop_complex;
                    ham=ham0+ham1;
                    hamm=kron(eye(n,n),ham)-kron(ham.',eye(n,n));
                    L=-1i*hamm-system.gamma;
                    sigss=system.gamma*sig_eq;  % steady state solutions for the denisty matrices
                    sigss=-L\sigss;
                    L=expm(L*wave.dt);
                    sigvec=sigss+L*(sigvec-sigss);
                end
                
                if ~isempty(options.detect) && experiment.detection(experiment.i)
                    sig=reshape(sigvec,n,n);
                    % computes expectation values
                    for j=1:ld
                        detected_signal(j,k+1)=sum(sum(options.detect{j}.*sig.'))/norms(j);
                    end
                    
                    % and stores density matrices if requested
                    if isfield(options,'store_density_matrix') && ~isempty(options.store_density_matrix) && options.store_density_matrix
                        sigmas{k+1}=sig;
                    end
                end
            end
            sig=reshape(sigvec,n,n);
        end
                
        % if phase cycling is requested, the density matrices and
        % expectation values are combined stored
        if size(wave.pcycle,1)>1
            pcycle_det=wave.pcycle(rr,2)/norm;
            loopstate=loopstate+pcycle_det*sig;
            if ~isempty(options.detect) && experiment.detection(experiment.i)
                if exist('pcycle_sig','var')
                    pcycle_sig=pcycle_sig+pcycle_det*detected_signal;
                else
                    pcycle_sig=zeros(size(detected_signal));                  
                    for lr=1:length(options.detect)
                        pcycle_sig=pcycle_det*detected_signal;
                    end
                end             
            end      
        end
        
        
    end
    
    % for phase cyclinge, the phase cycle results are rewritten to the
    % initial variables
    if size(wave.pcycle,1)>1
        state=loopstate;
        if ~isempty(options.detect) && experiment.detection(experiment.i)
            detected_signal = pcycle_sig;
        end
    else
        state = sig;
    end
       
else
    % Propagation in the case of a free evolution event
    
    % Setting up time axis for the two cases of direct and stepwise
    % evolution
    if ~experiment.stepwise_evolution(experiment.i)      
        t(2)=experiment.tp(experiment.i);
        dt=t(2);    
    else
        
        if ~isfield(experiment,'dt') || isempty(experiment.dt)  % if no timegrid is given, this creates a first guess from the length of the free evolution
            experiment.dt=experiment.tp(experiment.i);
        end
        
        % checks if the time grid is small enough to fullfill Nyquist
        if 1/experiment.dt<2*abs(system.highest_freq) 
            t=linspace(0,experiment.tp(experiment.i),round(experiment.tp(experiment.i)*2.5*system.highest_freq)+1); %If Nyquist is not met, a new time grid is created from the resonance frequency of the spinsystem
            dtn=t(2)-t(1);
            warning('timesteps too big for evolution, new timesteps have been created that fit the nyquist criterium.\n old: %6.0f ns; new:  %6.3f ns \n',experiment.dt,dtn);
            experiment.dt=dtn;
        else
            t=linspace(0,experiment.tp(experiment.i),ceil(experiment.tp(experiment.i)/experiment.dt)+1);
        end
        dt=experiment.dt;
    end
    
    sig=state; 
    ham=system.ham;
    
    
    if ~isempty(options.detect) && experiment.detection(experiment.i) 
        % setting up data structures for expectation values and density
        % matrices and computation of initial 
        
        ld=length(options.detect);
        detected_signal=zeros(ld,length(t));
        norms=zeros(1,length(options.detect));
        
        if isfield(options,'store_density_matrix') && ~isempty(options.store_density_matrix) && options.store_density_matrix
            sigmas=cell(1,length(t));
            sigmas{1}=sig;
        end
        
        for kk=1:length(options.detect),
            norms(kk) = sum(sum(options.detect{kk}.*options.detect{kk}));
            detected_signal(kk,1) = sum(sum(options.detect{kk}.*sig.'))/norms(kk);
        end
        
    end
    
    if options.relaxation == 0
        % propagation in Hilbert space
        U=expm(-1i*ham*dt);
        
        for k=2:length(t)
            
            sig=U*sig*U';
            
            
            if ~isempty(options.detect) && experiment.detection(experiment.i)
                % detection if requestesd
                for j=1:ld
                    detected_signal(j,k)=sum(sum(options.detect{j}.*sig.'))/norms(j);
                end
                % storing density matrices if requested
                if isfield(options,'store_density_matrix') && ~isempty(options.store_density_matrix) && options.store_density_matrix
                    sigmas{k}=sig;
                end
            end
        end
        
    else
        % setting up for propagation in Liouville space
        [n,~]=size(sig);
        sigvec=reshape(sig,n*n,1);
        sig_eq=reshape(system.eq,n*n,1);
        
        ham=kron(eye(n,n),ham)-kron(ham.',eye(n,n));
        L=-1i*ham-system.gamma;
        sigss=system.gamma*sig_eq;
        sigss=-L\sigss;
        L=expm(L*dt);
        
        for k=2:length(t)
            % propagation
            sigvec=sigss+L*(sigvec-sigss);
            sig=reshape(sigvec,n,n);
           
            if ~isempty(options.detect) && experiment.detection(experiment.i)
                % expectation values
                for j=1:ld
                    detected_signal(j,k)=sum(sum(options.detect{j}.*sig.'))/norms(j);
                end
                % density matrices
                if isfield(options,'store_density_matrix') && ~isempty(options.store_density_matrix) && options.store_density_matrix
                    sigmas{k}=sig;
                end
            end
        end
        
    end

    state=sig;
    tables=[];
    
end

% setting up eventually missing output
if isempty(options.detect) || ~experiment.detection(experiment.i)
    detected_signal=[];
end
if ~isfield(options,'store_density_matrix') || isempty(options.store_density_matrix) || ~options.store_density_matrix || ~experiment.detection(experiment.i)
    sigmas = []; 
end

% some indices, for processing of time axis and helpers of
% wave_propagation.m
experiment.t{experiment.i}=t+experiment.o;
experiment.o=experiment.o+t(end);
experiment.i=experiment.i+1;

end