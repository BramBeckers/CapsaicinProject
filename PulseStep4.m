%Pulse log processing script, step 4
%Author: Ali Gholamrezaei, University of Sydney


close all;
clear all;

%folder = ''; specify path to input files
files = dir(strcat(folder, '*.puls.txt'));

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
s_rmssd_deletion_base = [];

fdataname = strcat(folder, files(filenr).name);
casename = files(filenr).name;
string = strsplit(casename, '_');
case_n = regexp(string(3),'\d*','match');
condition = regexp(string(4),'\d*','match');

data = readtable(fdataname);
Rtimes = data.Rtimes;
IBIs = data.IBIs;
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

	if sum(error_final)>0
        IBIn = IBIs(error_final==0); % normal beats                    
        IBIs(error_final>0) = NaN;                    
        	if error_final(1)>0
            	IBIs(1)= IBIn(1);
            end                
            if error_final(end)>0
            	IBIs(end)= IBIn(end);
            end
        IBIc = fillmissing(IBIs,'pchip');
                
	elseif sum(error_final)==0
    	IBIc = IBIs;                
	end
                
%             time = (Rtimes(1):time_res/fs:Rtimes(end));
%             IBIint = interp1(Rtimes,IBIc,time','spline')';
% 
%             
%             y = filtfilt(but,a,IBIint);
%             
%             plot(Rtimes, IBIc)
%             hold on
%             plot(time, y);
%             
          
s_end = [s_end; floor((Rtimes(end)-(delay*time_res))/(w*time_res))];

for s = 1:s_end
    
    s_data = [];

    start_t = (delay+((s-1)*w))*time_res;
    end_t = (delay+(s*w))*time_res; 

    if s == 1
       start_t = start_t-(60*time_res);
    end
    
    s_beats = find(Rtimes>=start_t & Rtimes<=end_t);
    s_IBIs = IBIc(s_beats);
    
    s_errors = error_final(s_beats)'; 
    s_errors_per = round(sum(s_errors)*100/size(s_beats,1), 1);
         
%             sd = diff(s_IBIs);
%             
%             fault_sd = [];
%             
%             if s_errors_per>0
%             fault_sd = find(s_errors > 0);
%             
%             if fault_sd(1)>1
%             fault_sd = [fault_sd(1)-1;fault_sd];
%             end
%             sd(fault_sd) = 0;
%             nA = size(fault_sd, 1);
%             
%             else
%                 
%                 nA = 0;
%             
%             end
%                         
%             ssd = sd.^2;
%             ssd = sum(ssd);    
% 
%             N = size(s_IBIs, 1) - (nA+1);
%             mssd = ssd/N;
%             s_rmssd = sqrt(mssd);            
%             s_meanIBI = mean(s_IBIs(s_errors == 0));
            
            sd = diff(s_IBIs);
            ssd = sd.^2;
            mssd = mean(ssd);    
            s_rmssd_interpolation = round(sqrt(mssd),0);            
            s_meanIBI_interpolation = round(mean(s_IBIs),0); 
            
            
            time = (Rtimes(s_beats(1))/time_res:1/fs:Rtimes(s_beats(end))/time_res);
            
%             plot(Rtimes(s_beats(1):s_beats(end))/time_res, s_IBIs);
%             hold on
            
            s_IBIs = interp1(Rtimes(s_beats(1):s_beats(end))/time_res,s_IBIs,time','spline')';
            
%             plot(time, s_IBIs)
%             hold on
            
            s_IBIs = filtfilt(but,a,s_IBIs);
%             plot(time, s_IBIs)

%             fs = 4;
            
%             y = downsample(y, 1000/fs);

            window = min(300*fs, length(s_IBIs));
            noverlap = window/2;
            nfft = max(256,2^nextpow2(window));

    % (nfft/2)+1 if nfft is even, and (nfft+1)/2 if nfft is odd

            if rem(nfft, 2) == 0
                DFT = (nfft/2)+1;
            else
                DFT = (nfft+1)/2;
            end

            DFT = 0:0.001:2;

            order = 10;            
            [Pyy, f] = pburg(s_IBIs,order,DFT,fs, 'onesided');
            
%             peak = max(Pyy);
            
%             [Pyy, f] = pwelch(y,window,noverlap,DFT,fs,'PSD','onesided'); 


            VLF_ab = sum(Pyy(f >= 0.003 & f < 0.04));
            LF_ab = sum(Pyy(f >= 0.04 & f < 0.15));
            HF_ab = round(sum(Pyy(f >= 0.15 & f < 0.4)),0);
            
            Total_power = VLF_ab + LF_ab + HF_ab;
            
            VLF_p = VLF_ab*100/(Total_power); 
            LF_p = LF_ab*100/(Total_power);
            HF_p = round(HF_ab*100/Total_power,0);
            
            LF_nu = LF_ab*100/(LF_ab + HF_ab);
            HF_nu = round(HF_ab*100/(LF_ab + HF_ab),0);

            s_IBIs = IBIc(s_beats);
            s_IBIs(s_errors>0)=NaN;
            sd = diff(s_IBIs);
            ssd = sd.^2;
            mssd = mean(ssd,'omitnan');    
            s_rmssd_deletion = round(sqrt(mssd),0);            
            s_meanIBI_deletion = round(mean(s_IBIs(s_errors == 0),'omitnan'),0);            
    
%     s_data = [case_n condition ...
%         s start_t end_t ...
%         s_rmssd_interpolation s_rmssd_deletion ...
%         s_meanIBI_interpolation s_meanIBI_deletion ...
%         s_errors_per];
    if s==1
    s_rmssd_deletion_base = s_rmssd_deletion;
    end
    
    s_rmssd_deletion = s_rmssd_deletion-s_rmssd_deletion_base;

    
    s_data = [case_n condition ...
        s start_t end_t ...
        s_rmssd_deletion ...
        s_meanIBI_deletion ...
        HF_ab HF_p HF_nu ...
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

t = cell2table(allcasedata, 'VariableNames', {'case_n', 'condition', ...
    'sample_n', 'start_t', 'end_t', ...
    'rMSSD_deletion', ...
    'meanIBI_deletion', ...
    'HF_ab' 'HF_p' 'HF_nu', ...
    'sample_error_perc'});

writetable(t, 'D:\Bram\HRV_summary_5 min_adj.txt', 'delimiter', '\t');

% xlRange = 'E1';
% xlswrite('D:\Bram\HRV_summary_3 min.xls', allcasedata, xlRange);
