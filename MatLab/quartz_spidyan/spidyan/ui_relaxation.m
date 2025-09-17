function [system] = ui_relaxation( system, options )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    

    
    if system.spins==2
        nh=(system.sqn1*2+1)*(system.sqn2*2+1);
    elseif system.spins==1
        nh=(system.sqn1*2+1);
    end
    
    if ~isfield(system,'T1')
        system.T1=zeros(nh,nh);
    end
    
    [a,b]=size(system.T1);
    
    if a==b && a==1
        mat=ones(nh,nh);
        system.T1=system.T1*mat;
    elseif (a==b && a~=nh) || (a~=b)
        fprintf(1,'The input for T1 matrix was not recognized, using an empty matrix. \n');
        system.T1=zeros(nh,nh);
    end
    
    
    if ~isfield(system,'T2')
        system.T2=zeros(nh,nh);
    end
    [a,b]=size(system.T2);
    
    if a==b && a==1
        mat=ones(nh,nh);
        system.T2=system.T2*mat;
    elseif (a==b && a~=nh) || (a~=b)
        fprintf(1,'The input for T2 matrix was not recognized, using an empty matrix. \n');
        system.T2=zeros(nh,nh);
    end   
        
        

    
    f1 = figure('Position',[100 600 450 150],...
        'name','T1 table - close to continue','numbertitle','off');
    name =   {'1st state', '2nd state', '3rd state', 'etc...'};
    columnformat = cell(1,nh);
    [columnformat{:,:}] = deal('numeric');
    columneditable(1:nh) = true;
    
    callbackstring = 'system.T1=get(gco,''Data'');';
    Config = uitable('Data', system.T1,... 
                'ColumnName', name,...
                'ColumnFormat', columnformat,...
                'ColumnEditable', columneditable,...
                'celleditcallback', callbackstring,...
                'DeleteFcn',callbackstring,...
                'RowName',name);
    Config.Position(3) = Config.Extent(3);
    Config.Position(4) = Config.Extent(4);
    waitfor(f1);
    
    
    if isempty(system.T1)
        system.T1=zeros(nh,nh);
    end
 
    callbackstring = 'system.T2=get(gco,''Data'');';
    f2 = figure('Position',[100 600 450 150],...
        'name','T2 table - close to continue','numbertitle','off');
    Config2 = uitable('Data', system.T2,... 
                'ColumnName', name,...
                'ColumnFormat', columnformat,...
                'ColumnEditable', columneditable,...
                'celleditcallback',callbackstring,...
                'DeleteFcn',callbackstring,...
                'RowName',name);
    Config2.Position(3) = Config2.Extent(3);
    Config2.Position(4) = Config2.Extent(4);
    
    waitfor(f2);
    
    if isempty(system.T2)
        system.T2=zeros(nh,nh);
    end


end
