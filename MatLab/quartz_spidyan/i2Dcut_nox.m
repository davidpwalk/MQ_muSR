function ctlh = i2Dcut_nox(x1,x2,y2d,fignum)
% ctlh = i2Dcut(~x1,~x2,y2d,~fignum)
% creates an interactive cut through 2D data. x1 is the axis along the
% first dimension, x2 the axis along the second dimension. The first
% dimension will be used as x axis of the 1d cut. The figure number
% defaults to 1. Sometimes a call to i2Dcut_nox(y2d,fignum) makes live
% easier

% Andrin Doll, 2012-2015

% parse input

% figure number provided?
if ~exist('fignum','var')
    fignum = 1;
end

% x1,x2 omitted? (minimalistic useage)
if ~exist('y2d','var')
    % move first argument to y2d
    y2d = x1;
    % check for second argument
    if exist('x2','var')
        if length(x2) == 1
            fignum = x2;
        else
            error('For the minimalistic call, the second argument is the figure number')
        end
    end
    x1 = 1:size(y2d,1);
    x2 = 1:size(y2d,2);
end

figure(fignum)
clf
lh = line(x1,y2d(:,1)); % init with first variable
setappdata(gca,'line',lh);
for ii=1:length(x2)
    setappdata(gca,['dta' num2str(ii)],y2d(:,ii));
end
setappdata(gca,'x2',x2);
% ctlh.linedat = y2d;

% get handles and set figure props
fh = gcf;
ah = gca;

% This avoids flickering when updating the axis
set(fh,'doublebuffer','on');

% Generate constants for use in uicontrol initialization
pos=get(ah,'position');

% set(ah,'ylim',[-110 max(max(y2d))]);
% 
% % set grid from the actual value downto -120
% curr = get(ah,'YTick');
% set(ah,'Ytick',[-110:10:curr(end)]);
% grid on


% y text
ytext.pos = [pos(1)+pos(3)+0.005 pos(2) 0.06 0.05];
ytext.h = uicontrol(fh,'style','text', ...
    'units','normalized','position',ytext.pos,'String',num2str(x2(1)), ...
    'BackgroundColor',get(fh,'color'));
setappdata(gca,'ytext',ytext.h);

% yslide
yslide.pos=[pos(1)+pos(3)+0.02 pos(2)+0.1 0.03 pos(4)-0.2];
yslide.min = 1;
yslide.npts = 1000; % axis window dx by number of points
yslide.max = length(x2); % allows for overscaling
yslide.steps = length(x2)-1;
yslide.call=['set(getappdata(gca,''line''),''ydata'',getappdata(gca,[''dta'' num2str(floor(get(getappdata(gca,''yslide''),''value'')))]));'];
% now the x2 index...
yslide.call2=['set(getappdata(gca,''ytext''),''String'',num2str(subsref(getappdata(gca,''x2''),struct(''type'',''()'',''subs'',{{(floor(get(getappdata(gca,''yslide''),''value'')))}}))));' ];
% The big hack for matlab 2014b and newer: make a callback function out of
% that string (to make use of addlistener, instead of the old
% handle.listener)
eval(['yslide.call_fun = @(hObject,callbackdata) ' yslide.call ';']);
eval(['yslide.call_fun2 = @(hObject,callbackdata) ' yslide.call2 ';']);
yslide.h=uicontrol(fh,'style','slider',...
    'units','normalized','position',yslide.pos,...
    'callback',[yslide.call yslide.call2],'min',yslide.min,'max',yslide.max,'SliderStep',[1 1]./yslide.steps,'value',1);
ctlh.yslide = yslide;
% http://undocumentedmatlab.com/blog/continuous-slider-callback (in
% reversed order to avoid warnings..)
try  % R2014a and newer
    addlistener(yslide.h,'ContinuousValueChange',yslide.call_fun);
    addlistener(yslide.h,'ContinuousValueChange',yslide.call_fun2);
catch % R2013b and older
    addlistener(yslide.h,'ActionEvent',yslide.call_fun);
    addlistener(yslide.h,'ActionEvent',yslide.call_fun2);
end
setappdata(gca,'yslide',yslide.h);

set(gcf,'Toolbar','figure')
% keyboard
