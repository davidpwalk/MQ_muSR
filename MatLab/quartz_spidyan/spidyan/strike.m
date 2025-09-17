function [signal_dc,filter] = strike(signal_unfiltered, t, down, options)
% [signal_filtered,filter] = strike(signal_unfiltered, t, down, options)
% 
% removes oscillation from a given detected signal, Sz, S1z and S2z are not
% filtered by default, since only signals with x or y component contain
% oscillations. If down conversion is not wanted for a certain signal (eg. if the
% resonance frequency of the second spin is much smaller) the detection
% operator as string or the number of the signal can be specified. Used to 
% transfer signal from rotating frame to lab frame. options.det_op must be
% present for strike to function.
% 
% Input:
% signal_unfiltered detected signal as it is returned by wave_propagation
% t                 time axis of the detected signal
% down              frequency for downconversion in GHz
%          
%                   
% .options      structure containing additional parameters
% .det_op       cell containing the detection operators as string (S1z,S2z,
%               Sz,...) or in matrix matrix form. For more details see
%               documentation, can be empty.
% .no_dc        cell which containing the signals not to be filtered as
%               string or index of detection operator, optional.
% .cutoff_freq  cutoff frequency of the applied filter for down conversion
%               in MHz. Usually 200 MHz is a good value to start with. 
% .LO           frequency of the local oscillator [GHz], required for
%               computations in the lab frame.
% .labframe     0 or 1, default is 0. If true SPIDYAN simulates the
%               experiment in a local frame rotating at
%               options.LO+system.nu01.
% .display_filter   0 or 1, optional, default is false. If true displays the
%                   frequency response of the filter used for downconversion.
% 
% Output:              
% signal_filtered       returns the (complex) detected signal without 
%                       oscillations.
% filter         structure with the filter.
% 
% St. Pribitzer, 2014

if ~exist('options','var') 
    options=[];
end

cutoff=200; % default cutoff frequency is 200 MHz

if isfield(options,'cutoff_freq')
     cutoff=options.cutoff_freq;
end

if isfield(options,'det_op')
     detect_op=options.det_op;
end


kk=1;
for k=1:length(detect_op)
    % checks what elements of the signals are to be filtered and creates
    % a vector v which contains the index of the signals
    a=1;
    
    if ischar(detect_op{k}) % signals of which the detection operator 
        % contains x,y,m or p elements are filter, unless excluded with
        % options.no_dc
        
        if ~any(detect_op{k}=='x') && ~any(detect_op{k}=='y') ...
                && ~any(detect_op{k}=='p') && ~any(detect_op{k}=='m')
            
            a=0;
            
        else
            
            if isfield(options,'no_dc') && ~isempty(options.no_dc)
                for r=1:length(options.no_dc)
                    if ischar(options.no_dc{r})
                        if strcmp(detect_op{k},options.no_dc{r})
                            a=0;
                        end
                    elseif options.no_dc{r}==k
                        a=0;
                    end
                end       
            end 
            
        end
        
        if a==1
            v(kk)=k; %#ok<AGROW>
            kk=kk+1;
        end 
        
    else
        % if the detection operator is given as matrix, SPIDYAN filters
        % those signals with off-diagonal elements in the detection
        % operator, unless excluded with options.no_dc
        
        [n, ~]=size(detect_op{k}); 
        c=0;
        
        for mm=1:n
            for nn=1:n
                if mm~=nn && detect_op{k}(mm,nn)~=0
                    c=c+1;
                end
            end
        end
        
        if c~=0
            v(kk)=k; %#ok<AGROW>
            kk=kk+1;
        end 
        
    end
end

      
if exist('v','var')
    
    signal=zeros(length(signal_unfiltered),length(v));
    
    % signals which qualified for down conversion are stored in a new matrix
    for k=1:length(v)
        signal(:,k)=signal_unfiltered(v(k),:);
    end
    
    dt=t(2)-t(1);
    
    % if simulation was carried out in the labframe, the frequency to be
    % filtered is adapted to it
    if isfield(options,'labframe') && ~isempty(options.labframe) && options.labframe && isfield(options,'LO') && ~isempty(options.LO)
        down=down+options.LO;
    end
    
    % oscillation to removed
    rf = exp(-1i*(2*pi*(t-dt/2)*down));

    % mulitplication of unfiltered signals with the oscillation
    [x,y]=size(signal);
    signal2 = bsxfun(@times,signal,rf.'); 
    cutoff_freq=cutoff;
    
    % butterworth filter in timedomain signal -> sig3
    cutoff=cutoff/(options.awg.s_rate/2)/1000;
    
    % filter is a lowpass butterworth filter of 3rd order
    N=3;
    [B,A] = butter(N,cutoff,'low');
    
    % creates frequency responce of filter
	[H,F]=freqz(B,A,512,12);
    
    filter.H=H;
    filter.F=F;
    
    if isfield(options,'display_filter') && options.display_filter
        % if true, frequecy respond of the filter is displayed
        figure(99);clf; 
        plot(F*1000,abs(H).^2) % signal is filtered twice with filtfilt, therefore ^2
        xlim([0 5*cutoff_freq])
        xlabel('f [MHz]')
        ylabel('H [a.u.]')
        title('Frequency response of filter')     
    end

    % signal has to be expanded by the response time of the filter at the
    % beginning and the end to minimize artifacts
    taufilter=2*cutoff*N;
    zerofill=5*taufilter;
    nzerofill= round(zerofill/dt);
    zerosig=zeros(x+2*nzerofill,y);
    
    zerosig(nzerofill+1:nzerofill+length(signal2),:)=signal2(:,:);
    
    % signal is filtered with filtfilt, which provides almost zero phase
    % distortion
    zerofilt = filtfilt(B,A,zerosig);
    
    sig3=zerofilt(nzerofill+1:nzerofill+length(signal2),:);
    sig3=2*sig3;

    [h, ~]=size(signal_unfiltered);

    signal_dc=zeros(size(signal_unfiltered));
    j=1;
    
    for k=1:h
        % restoration of all signals into one matrix again
        if j <= length(v)
           if all(signal_unfiltered(k,:)==0) % checks if the signal is only
                    % zero, if yes, the unfiltered signal is returned
               signal_dc(k,:)=signal_unfiltered(k,:);
           elseif v(j)==k
               signal_dc(k,:)=sig3(:,j);
               j=j+1;
           else
               signal_dc(k,:)=signal_unfiltered(k,:);
           end
        else
            signal_dc(k,:)=signal_unfiltered(k,:);
        end
    end
    
else
    signal_dc=signal_unfiltered;
    filter=[];
end

end



