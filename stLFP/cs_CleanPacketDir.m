clear; close all;

packetPath = 'F:\S\Packets\';
folders = dir(packetPath); folders = folders(3:end);

for f = 1:length(folders)
   filepath = [packetPath,folders(f).name,'\'];
   cd(filepath)
   delete *.ps;
   delete *.log;
end