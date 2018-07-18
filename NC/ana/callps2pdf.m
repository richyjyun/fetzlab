function callps2pdf(packet,replace)

if(~exist('replace','var'))
    replace = 1;
end

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