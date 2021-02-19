%Resp log processing script, step 2
%Author: Ali Gholamrezaei, University of Sydney


close all;
clear all;

%folder = ''; specify path to input files
files = dir(strcat(folder, '*.resp'));

allcasedata = [];

for filenr = 1:size(files,1)

fname = strcat(folder, files(filenr).name);

fid=fopen(fname);

alldata = textscan(fid,'%s');

for n = 1:size(alldata{1},1)
   if strcmp(alldata{1}(n),'RESP_SAMPLES_PER_SECOND')
       RESP_SAMPLES_PER_SECOND = alldata{1}(n+2);
   end
   if strcmp(alldata{1}(n),'RESP_SAMPLE_INTERVAL')
       RESP_SAMPLE_INTERVAL = alldata{1}(n+2);
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
   if strcmp(alldata{1}(n),'RESP_SAMPLE_INTERVAL')
       startsignal = n+3;
   end
   if strcmp(alldata{1}(n),'ACQ')
       endsignal = n-1;
   end
end

RESP_SAMPLES_PER_SECOND = split(RESP_SAMPLES_PER_SECOND, ';');
RESP_SAMPLES_PER_SECOND = str2double(RESP_SAMPLES_PER_SECOND{1,1});
RESP_SAMPLE_INTERVAL = str2double(RESP_SAMPLE_INTERVAL);
LogStartMDHTime = str2double(LogStartMDHTime);
LogStopMDHTime = str2double(LogStopMDHTime);
LogStartMPCUTime = str2double(LogStartMPCUTime);
LogStopMPCUTime = str2double(LogStopMPCUTime);

fdataname = strcat(fname, '.txt');
if isfile(fdataname)
data = readtable(fdataname);

errors = data.errors;
error_final = sum(errors>=1)*100/size(data,1);
error_final = round(error_final,3);

rows = errors>=1;

EXPtimes = data.EXPtimes;
INSPtimes = data.INSPtimes;
EXPtimes(rows)=NaN;
INSPtimes(rows)=NaN;
RESPtime_EXP = round(nanmean(diff(EXPtimes)),0);
RESPtime_INSP = round(nanmean(diff(INSPtimes)),0);
RESPtime_Count_Hz = size(errors(errors~=1),1)/(EXPtimes(end)/1000);
RESPtime_Count_t = 1000/RESPtime_Count_Hz;

casename = files(filenr).name;
casedata = [];
casedata = {casename, RESP_SAMPLES_PER_SECOND, RESP_SAMPLE_INTERVAL, LogStartMDHTime, ...
   LogStopMDHTime, LogStartMPCUTime, LogStopMPCUTime, RESPtime_EXP, RESPtime_INSP, RESPtime_Count_Hz, RESPtime_Count_t, error_final};

allcasedata = [allcasedata; casedata];
end
end

t = cell2table(allcasedata, 'VariableNames', {'casename' 'RESP_SAMPLES_PER_SECOND' 'RESP_SAMPLE_INTERVAL' ...
    'LogStartMDHTime' 'LogStopMDHTime' 'LogStartMPCUTime' 'LogStopMPCUTime' 'RESPtime_EXP' ...
    'RESPtime_INSP' 'RESPtime_Count_Hz' 'RESPtime_Count_t' 'error_final'});

writetable(t, 'D:\Bram\summary_resp.txt', 'delimiter', '\t');

xlRange = 'E1';
xlswrite('D:\Bram\summary_resp.xls', allcasedata, xlRange);