%Pulse log processing script, step 2
%Author: Ali Gholamrezaei, University of Sydney


close all;
clear all;

%folder = ''; specify path to input files
files = dir(strcat(folder, '*.txt'));

cases = [];

for filenr = 1:size(files,1)

fname = strcat(folder, files(filenr).name);

data = readtable(fname);

            Rtimes = data.Rtimes;

            if Rtimes(1)<200
                data=data(2:end,:);
                writetable(data, fname, 'Delimiter', '\t');
                
                cases = [cases; fname];
                
            end

end