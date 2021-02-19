%Pulse log processing script, step 3
%Author: Ali Gholamrezaei, University of Sydney


close all;
clear all;

%folder = ''; Specify path to input files
files = dir(strcat(folder, '*.puls'));
allcasedata = [];
for filenr = 1:size(files,1)

fname = strcat(folder, files(filenr).name);

fid=fopen(fname);

alldata = textscan(fid,'%s');

for n = 1:size(alldata{1},1)
   if strcmp(alldata{1}(n),'PULS_SAMPLES_PER_SECOND')
       PULS_SAMPLES_PER_SECOND = alldata{1}(n+2);
   end
   if strcmp(alldata{1}(n),'PULS_SAMPLE_INTERVAL')
       PULS_SAMPLE_INTERVAL = alldata{1}(n+2);
   end
   if strcmp(alldata{1}(n),'LogStartMDHTime:')
       LogStartMDHTime = alldata{1}(n+1);
   end
   if strcmp(alldata{1}(n),'LogStopMDHTime:')
       LogStopMDHTime = alldata{1}(n+1);
   end
   if strcmp(alldata{1}(n),'LogStartMPCUTime:')
       LogStartMPCUTime = alldata{1}(n+1);
   end
   if strcmp(alldata{1}(n),'LogStopMPCUTime:')
       LogStopMPCUTime = alldata{1}(n+1);
   end   
   if strcmp(alldata{1}(n),'PULS_SAMPLE_INTERVAL')
       startsignal = n+3;
   end
   if strcmp(alldata{1}(n),'ACQ')
       endsignal = n-1;
   end
end

PULS_SAMPLES_PER_SECOND = split(PULS_SAMPLES_PER_SECOND, ';');
PULS_SAMPLES_PER_SECOND = str2double(PULS_SAMPLES_PER_SECOND{1,1});
PULS_SAMPLE_INTERVAL = str2double(PULS_SAMPLE_INTERVAL);
LogStartMDHTime = str2double(LogStartMDHTime);
LogStopMDHTime = str2double(LogStopMDHTime);
LogStartMPCUTime = str2double(LogStartMPCUTime);
LogStopMPCUTime = str2double(LogStopMPCUTime);
    
    
fdataname = strcat(fname, '.txt');

data = readtable(fdataname);

sum_criteria_th = data.sum_criteria_th;
errors = data.errors;
error_final = [];
for beat = 1:size(data,1)
      if sum_criteria_th(beat)==1 | errors(beat)==1
          error_final(beat) = 1;
      else 
          error_final(beat) = 0;
          
      end

end

error_final = sum(error_final)*100/size(data,1);
error_final = round(error_final,3);

casename = files(filenr).name;
casedata = [];
casedata = {casename, PULS_SAMPLES_PER_SECOND, PULS_SAMPLE_INTERVAL, LogStartMDHTime, ...
   LogStopMDHTime, LogStartMPCUTime, LogStopMPCUTime, error_final};

allcasedata = [allcasedata; casedata];

end

t = cell2table(allcasedata, 'VariableNames', {'casename', 'PULS_SAMPLES_PER_SECOND', ...
    'PULS_SAMPLE_INTERVAL', 'LogStartMDHTime', ...
   'LogStopMDHTime', 'LogStartMPCUTime', 'LogStopMPCUTime', 'error_final'});

writetable(t, 'D:\Bram\summary.txt', 'delimiter', '\t');

xlRange = 'E1';
xlswrite('D:\Bram\summary.xls', allcasedata, xlRange);
