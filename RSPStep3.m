%Resp log processing script, step 3
%Author: Ali Gholamrezaei, University of Sydney


close all;
clear all;

%folder = ''; specify path to input files
files = dir(strcat(folder, '*.resp.txt'));

delay = 25+60; %initial delay in sec 
w = 5*60; %windows duration in sec
time_res = 1000; %time resolution in data file, 1000 = msec

fs = 4; %sampling freq for spectral analysis
N = 1;  %first order
cutoff_Hz = 0.015;  %should be 3dB down at the cutoff
[but,a]=butter(N,cutoff_Hz/(fs/2),'high'); 
            
allcasedata = [];

s_end = [];

for filenr = 1:size(files,1)

casedata = [];
    
fdataname = strcat(folder, files(filenr).name);
casename = files(filenr).name;
string = strsplit(casename, '_');
case_n = regexp(string(3),'\d*','match');
condition = regexp(string(4),'\d*','match');

data = readtable(fdataname);
EXPtimes = data.EXPtimes;
RESPampl = data.RESPampl;
EXPIntervals = data.EXPIntervals;
INSPIntervals = data.INSPIntervals;
errors = data.errors;

s_end = [s_end; floor((EXPtimes(end)-(delay*time_res))/(w*time_res))];

for s = 1:s_end
    
    s_data = [];

    start_t = (delay+((s-1)*w))*time_res;
    end_t = (delay+(s*w))*time_res; 

    if s == 1
       start_t = start_t-(60*time_res);
    end
    
    s_RSP = find(EXPtimes>=start_t & EXPtimes<=end_t);
    s_errors = errors(s_RSP);
    s_errors_per = round(sum(s_errors)*100/size(s_RSP,1), 1);

    s_RSP = find(EXPtimes>=start_t & EXPtimes<=end_t & errors==0);

    s_RESPampl = round(mean(RESPampl(s_RSP)),1);
    s_EXPIntervals = round(mean(EXPIntervals(s_RSP)),1);
    s_INSPIntervals = round(mean(INSPIntervals(s_RSP)),1);
         
    s_BR = round(size(s_RSP,1) / (((end_t-start_t)/time_res)/60),1);
    

    s_data = [case_n condition ...
        s start_t end_t ...
        s_BR s_EXPIntervals s_INSPIntervals s_RESPampl...
        s_errors_per];
    
    casedata = [casedata; s_data];

end

allcasedata = [allcasedata; casedata];

end

% t = cell2table(allcasedata, 'VariableNames', {'case_n', 'condition', ...
%     'sample_n', 'start_t', 'end_t', ...
%     'rMSSD_interpolation', 'rMSSD_deletion', ...
%     'meanIBI_interpolation', 'meanIBI_deletion', ...
%     'sample_error_perc'});

t = cell2table(allcasedata, 'VariableNames', {'case_n' 'condition' ...
        's' 'start_t' 'end_t' ...
        's_BR' 's_EXPIntervals' 's_INSPIntervals' 's_RESPampl' ...
        's_errors_per'});

writetable(t, 'D:\Bram\RSP_summary_5 min.txt', 'delimiter', '\t');

% xlRange = 'E1';
% xlswrite('D:\Bram\HRV_summary_3 min.xls', allcasedata, xlRange);
