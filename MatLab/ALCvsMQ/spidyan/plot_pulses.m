function plot_pulses(experiment,options)
% plot_pulses(experiment,options)
% 
% This function uses the experiment structure from SPIDYAN to plot the
% pulses. Several configurable parameters allow the user to select the
% pulses which are to be plotted, whether to plot in the time or frequency
% domain or both and if he wants the imaginary part to be displayed or not.
% 
% Parameters of options
% .plot_domain          'time', 'frequency' or 'both'. Selects what plots
%                              to display, default is to 'time'
% .plot_imaginary_part  0/1, plot imaginary part of wave, default is 0
% .plot_pulse           boolean vector, ordering corresponds to
%                       experiment.tp, 1 for plotting, 0 for no plotting,
%                       if empty or not defined, all pulses are plotted.
% 
% St. Pribitzer, 2015

% Setting up of default parameters if necessary
if ~isfield(options,'plot_domain') || isempty(options.plot_domain)
    options.plot_domain = 'time';
end
if ~isfield(options,'plot_imaginary_part') || isempty(options.plot_imaginary_part)
    options.plot_imaginary_part = 0;
end

% indices
r = 1;
toplot = 0;

% processing of boolean vector, not pretty, could be combined with the
% plotting loop
if isfield(options,'plot_pulses') && ~isempty(options.plot_pulses)
    for k = 1: length(options.plot_pulses)
        if k <= length(experiment.pulse) && options.plot_pulses(k) && ~isempty(experiment.pulse{k})
            toplot(r)=k;
            r = r+1;
        end
    end         
elseif ~isfield(options,'plot_pulses') || isempty(options.plot_pulses)
    for k = 1: length(experiment.pulse)
        if ~isempty(experiment.pulse{k})
            toplot(r)=k;
            r  = r+1;
        end
    end
end

% This is were the plotting happens by calling either of the two local
% functions plotfrequency and plottime
if length(toplot)>=1 && toplot(1) ~= 0
    for k = 1: length(toplot)
        switch options.plot_domain
            case 'time'
                plottime(experiment.pulse{toplot(k)},options)
            case 'frequency'
                plotfrequency(experiment.pulse{toplot(k)},options)
            case 'both'
                plottime(experiment.pulse{toplot(k)},options)
                plotfrequency(experiment.pulse{toplot(k)},options)
        end
    end
else
    warning('No pulses selected to plot')
end

end

function plotfrequency(pulse,options)
% local function, which plots the pulse in the frequency domain, requires
% the funtion dft, provided with SPIDYAN
opt.force_complete = 1;
[nu,spc] = dft(pulse.t,pulse.binary-mean(pulse.binary),opt);
spc = spc/max(spc);

figure


if options.plot_imaginary_part
    hold on
    plot(nu,real(spc))
    plot(nu,imag(spc))
else
    plot(nu,abs(spc))
end
xlabel('\nu [GHz]')
ylabel('Int. [a.u.]')
end


function plottime(pulse,options)
% local function to that plots the pulse in the time domain
figure
plot(pulse.t,real(pulse.binary))

if options.plot_imaginary_part
    hold on
    plot(pulse.t,imag(pulse.binary))
end
xlabel('t [ns]')
ylabel('Output level')

end