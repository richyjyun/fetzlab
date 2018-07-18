function callps2pdf(packet,replace)
% Summary
%   Calls ps2pdf function for converting ps files to pdf. 
%
% Inputs
%   packet - packet name
%   replace - whether or not to replace the current pdf
%
% RJY 06/22/2018

    % Set variables
    if(~exist('replace','var'))
        replace = 1;
    end

    % Error checking
    if(~exist(packet,'file'))
        error('PS FILE DOES NOT EXIST');
    end

    if(exist([packet(1:end-1),'df'],'file'))
        if(replace)
            delete([packet(1:end-1),'df']);
        else
            error('FILE ALREADY EXISTS');
        end
    end

    % Call ps2pdf
    try
        ps2pdf('psfile',packet,'pdffile',[packet(1:end-1),'df'],...
            'gscommand','C:\Program Files (x86)\gs\gs9.22\bin\gswin32.exe',...
            'gsfontpath','C:\Program Files (x86)\gs\gs9.22\lib',...
            'gslibpath','C:\Program Files (x86)\gs\gs9.22\lib');
        delete(packet);
    catch
        error('CONVERSION TO PDF FAILED');
    end

end