%Resp log processing script, step 1
%Author: Ali Gholamrezaei, University of Sydney


close all;
clear all;

%folder = ''; specify path to input files
files = dir(strcat(folder, '*.resp'));

Fs = 1000; 
Order = 1;  %first order
Fc = [0.1 1];  %should be 3dB down at the cutoff
[but,a] = butter(Order,Fc/(Fs/2),'bandpass'); 

MinRESPePD = 2000;
MinRESPePH = 90; %max(RESPeraw)/div;

%errors:
minVP = 500;
maxVP = 10000;
minEXPdu = 500;
maxRESPIntdiff = 1000;
maxRESPInterval = 8000;
maxpksdiff = -90;
maxampldiff = -90;

forward = 1000;
backward = 1000;

sum_criteria_t = 5;

movW = 15;
movT = 350;

for filenr = 36%:size(files,1)

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

signal = alldata{1}(startsignal:endsignal);
signal = str2double(signal);

time = 1:1000/RESP_SAMPLES_PER_SECOND:length(signal)*1000/RESP_SAMPLES_PER_SECOND;
time = time';

% plot(time, signal);
% hold on

system_triggers = find(signal(:,1)>=5000);
RESPtimes_system = 1000/RESP_SAMPLES_PER_SECOND * system_triggers;
RESPIntervals_system = [NaN;diff(RESPtimes_system)];

signal_edited = signal;

signal_edited(system_triggers) = NaN;

signal_edited = fillmissing(signal_edited, 'spline');  % or pchip

% signal_edited = [NaN; diff(signal_edited)];
% plot(time, signal_edited); 

time2 = 1:1000/Fs:time(end);
signal_int = interp1(time,signal_edited,time2,'spline')'; %

signal_filtered = filtfilt(but,a,signal_int);

yyaxis left
plot(time2, signal_int);
yyaxis right
plot(time2, signal_filtered);
refline(0,0);

RESPeraw = signal_filtered; %signal_filtered;

            [pks, RESPtimes] = findpeaks(RESPeraw, 'MinPeakHeight', MinRESPePH, 'MinPeakDistance', MinRESPePD);
            
            if RESPtimes(1)<200
                RESPtimes=RESPtimes(2:end);
                pks=pks(2:end);
            end
            
            RESPIntervals = [NaN; diff(RESPtimes)];
%             pksdiff = [NaN; abs(diff(pks))];
%             pksdiff = pksdiff/nanmean(pksdiff);
            
            SD = [NaN; diff(RESPIntervals)];
            
            SDab = abs(SD);
            
            SDm = movmean(SDab, movW);
            SDmdiff = 100*abs(SDm-SDab)./SDm;

            V = [NaN];
            I = [NaN];
            VP = [NaN];

            
            
            for   RESPNr = 1:size(RESPtimes,1)-1
                v = [];
                i = [];
                [v, i] = min(RESPeraw(RESPtimes(RESPNr):RESPtimes(RESPNr+1)));
                i = RESPtimes(RESPNr) + i;
                vp = RESPtimes(RESPNr+1) - i;
                
                V = [V; v];
                I = [I; i];
                VP = [VP; vp];
            end
            
            INSPIntervals = [NaN; diff(I)];

            
            ampl = pks-V;
            ampldiff = [];
            pksdiff = [];
            RESPIntdiff = [];
            RESPIntervals3 = [];
            RESPtimes2 = RESPtimes;
            RESPIntervals2 = RESPIntervals;
            for RESPNr = 1: size(RESPtimes,1)
%             
% %                 RESPtimes2 = Rtimes;
                RESPtimes2(RESPNr)=NaN;
                RESPtimes2 = fillmissing(RESPtimes2, 'spline');  % or pchip
                RESPIntdiff(RESPNr) = abs(RESPtimes(RESPNr) - RESPtimes2(RESPNr));
%                 
                if RESPNr > 2 & RESPNr < size(RESPtimes,1)
                RESPIntervals2(RESPNr)=NaN;
                RESPIntervals2 = fillmissing(RESPIntervals2, 'spline');  % or pchip
                RESPIntdiff(RESPNr) = abs(RESPIntervals(RESPNr) - RESPIntervals2(RESPNr));
%                 
                RESPIntervals3(RESPNr) = abs(RESPIntervals(RESPNr) - mean(RESPIntervals(RESPNr-1),RESPIntervals(RESPNr+1)));
                
                pksdiff(RESPNr) = (pks(RESPNr)-(pks(RESPNr-1)+pks(RESPNr+1)))*100/...
                    (pks(RESPNr-1)+pks(RESPNr+1));
                ampldiff(RESPNr) = (ampl(RESPNr)-(ampl(RESPNr-1)+ampl(RESPNr+1)))*100/...
                    (ampl(RESPNr-1)+ampl(RESPNr+1));
                   
                else
                    RESPIntervals2(RESPNr)=NaN;
                    RESPIntervals3(RESPNr)=NaN;
                    pksdiff(RESPNr)=NaN;
                    ampldiff(RESPNr)=NaN;
                end                

            end

            RESPIntervals3=RESPIntervals3';
%             RESPIntervals3 = [NaN;NaN;RESPIntervals3];
            RESPIntdiff = RESPIntdiff';
%             RESPIntdiff = [NaN;RESPIntdiff];
            pksdiff = pksdiff';
            ampl = ampl';
%             sum_criteria =[];
            
%             RIntdiffsort = sort(RIntdiff, 'descend');
%             percentile5 = prctile(RIntdiff,[1:100], 'all'); 
            
%             L = log(RIntdiff);
            
            errors = []; 

            for   RESPNr = 1:size(RESPtimes,1)-1
            
                                           
            if (RESPIntdiff(RESPNr)>maxRESPIntdiff & VP(RESPNr)>maxVP) ...
                   | VP(RESPNr)<minVP | (VP(RESPNr)>maxVP & RESPIntdiff>maxRESPIntdiff) ...
                   | I(RESPNr+1)-RESPtimes(RESPNr)<minEXPdu ... % | V(beatNr)<-877 % | SDmdiff(beatNr+1)>movT 
                   | pksdiff(RESPNr) < maxpksdiff ...
                   | ampldiff(RESPNr) < maxampldiff ...
                   
                errors(RESPNr)=1;
            else 
                errors(RESPNr)=0;
                
            end
                        
            flatf = [NaN];
            flatb = [NaN];
                
                if RESPtimes(RESPNr) > forward
                flatb = signal(((RESPtimes(RESPNr)-backward)*RESP_SAMPLES_PER_SECOND/1000):...
                    (RESPtimes(RESPNr)*RESP_SAMPLES_PER_SECOND/1000));
                end
                
                flatf = signal((RESPtimes(RESPNr)*RESP_SAMPLES_PER_SECOND/1000):...
                    ((RESPtimes(RESPNr)+forward)*RESP_SAMPLES_PER_SECOND/1000));
                 
                if isnan(flatb)==0 & sum(diff(flatb))==0
                errors(RESPNr)=3;
                end
                
                if isnan(flatf)==0 & sum(diff(flatf))==0 
                errors(RESPNr)=3;
                end 
                 
            
%             sum_criteria(beatNr) = RIntdiff(beatNr)>maxRIntdiff + VP(beatNr)>maxVP + ... 
%                 VP(beatNr)<minVP + IBIs3(beatNr)>maxIBIs3 + SDmdiff(beatNr)>movT;
            
            end

            errors = [errors'; NaN];
            
%             sum_criteria =sum_criteria';
%             sum_criteria = [NaN;sum_criteria];
            
%             SDmdiff2 = [SDmdiff(2:end);NaN];

%             sum_criteria = ((RESPIntdiff>maxRIntdiff)*1) + (VP>maxVP) + (VP>400) + (VP>700) + ...
%                 ((RESPIntdiff>maxIBIIntdiff)*1) + (RESPIntdiff>250) + ... 
%                 (VP<minVP) + ((RESPIntervals3>maxIBIs3)*1) + (RESPIntervals3>300) + (RESPIntervals3>500) + ...
%                 (SDmdiff2>movT);
            
            Time = 1:size(RESPeraw,1);
                Time = downsample(Time, 5);
                RESPeraw = downsample(RESPeraw, 5);
              
                nf = 1;
%%
            figure(nf);
            ax1 = subplot(4,1,1);
            yyaxis left
            plot(Time, RESPeraw);
            hold on
            plot(RESPtimes, pks, 'o', 'Color', 'b');
            hold on
            plot(I, V, 'o', 'Color', 'r');
            
            if nanmean(errors)~=0
            hold on
            line([RESPtimes(errors>=1,1) RESPtimes(errors>=1,1)]', ylim, 'Marker','none','LineStyle','-','Color',[1 0 0]);
            end
%             if nanmean(sum_criteria)~=0 & isempty(RESPtimes(sum_criteria>=sum_criteria_t,1))~=1
% %             hold on
%             line([RESPtimes(sum_criteria>=sum_criteria_t,1) RESPtimes(sum_criteria>=sum_criteria_t,1)]', ylim, ... 
%                 'Marker','none','LineStyle','--','Color', 'K');
%             end
            ylim([-3000 3000]);
            xlim([min(Time), max(Time)]);
            title('RESPe');

            yyaxis right
            plot(time, signal, 'LineStyle', ':');
            
            ax2 = subplot(4,1,2);
            yyaxis left 
            plot(RESPtimes, RESPIntervals, '-o');
            hold on
            plot(RESPtimes, RESPIntervals2, '-*');
            hold on
            plot(RESPtimes_system, RESPIntervals_system, '-+', 'Color', 'K');
            
            
            ref1500 = refline(0 , 10000);
            ref500 = refline(0 , 2000);
            ref3std = refline(0 , (nanmean(RESPIntervals) + 3*nanstd(RESPIntervals)));
            ref3std.Color = 'b';
            ref3_std = refline(0 , (nanmean(RESPIntervals) - 3*nanstd(RESPIntervals)));
            ref3_std.Color = 'b';
            if nanmean(errors)~=0
            hold on
            line([RESPtimes(errors>=1,1) RESPtimes(errors>=1,1)]', ylim, 'Marker','none','LineStyle','-','Color',[1 0 0]);
            end
%             if nanmean(sum_criteria)~=0 & isempty(RESPtimes(sum_criteria>=sum_criteria_t,1))~=1
% %             hold on
%             line([RESPtimes(sum_criteria>=sum_criteria_t,1) RESPtimes(sum_criteria>=sum_criteria_t,1)]', ylim, ... 
%                 'Marker','none','LineStyle','--','Color', 'K');
%             end                       
            ylim([min(RESPIntervals)-2000, max(RESPIntervals)+2000]);
            xlim([min(Time), max(Time)]);
%             title('IBI');
            ylabel('msec');
                        
            yyaxis right
%             plot(RESPtimes, SD, '*');
            hold on
            plot(RESPtimes, SDmdiff, '-o');
%             refline(0 , movT);
            ylim([min(SDmdiff)-1000, max(SDmdiff)+1000]);
            xlim([min(Time), max(Time)]);
            title('IBI, SD, SD diff from movemean');
            ylabel('msec / % from movmean');
     
            ax3 = subplot(4,1,3);
            yyaxis left
            plot(RESPtimes, VP, 'o');
            refline(0 , 2000);
            ylabel('msec');
            if nanmean(errors)~=0
            hold on
            line([RESPtimes(errors>=1,1) RESPtimes(errors>=1,1)]', ylim, 'Marker','none','LineStyle','-','Color',[1 0 0]);
%             
%             
            end
%             if nanmean(sum_criteria)~=0 & isempty(RESPtimes(sum_criteria>=sum_criteria_t,1))~=1
% %             hold on
%             line([RESPtimes(sum_criteria>=sum_criteria_t,1) RESPtimes(sum_criteria>=sum_criteria_t,1)]', ylim, ... 
%                 'Marker','none','LineStyle','--','Color', 'K');
%             
%             
%             end
            yyaxis right
            plot(RESPtimes, RESPIntdiff, '^');
            hold on
            plot(RESPtimes, RESPIntervals3, '+', 'Color', 'K');
%             
%             
            refline(0 , 1000);
            refline(0, (nanmean(RESPIntdiff)+(3*nanstd(RESPIntdiff))));
%             title('D2S time, IBI diff if interpolated');
%             ylabel('msec');
            
            ax4 = subplot(4,1,4);
%             yyaxis left
            plot(RESPtimes, pksdiff, '-o');
            if nanmean(errors)~=0
            hold on
            line([RESPtimes(errors>=1,1) RESPtimes(errors>=1,1)]', ylim, 'Marker','none','LineStyle','-','Color',[1 0 0]);
            end
            
%             yyaxis right
            plot(RESPtimes, ampldiff, '-o')
            refline(0 , -90);
% %             ylabel('msec');

            linkaxes([ax1,ax2,ax3, ax4], 'x');
%%

end

	data = [RESPtimes pks I V VP RESPIntervals INSPIntervals ...
        RESPIntdiff RESPIntervals3 errors];
        
        data = round(data, 1);
        
    data = array2table(data, 'variablenames', {'EXPtimes' 'EXPampl' ...
        'INSPtimes' 'INSPampl' 'RESPampl' 'EXPIntervals' 'INSPIntervals' ...
        'Interpoldiff' 'movediff' 'errors'});
    
    writetable(data, strcat(fname, '.txt'), 'Delimiter', '\t');
    
    saveas(nf, strcat(fname, '_.fig'));
    
    rows = find(errors>=1);
    
    data(rows, :)
    
    (sum(errors>=1))*100/size(RESPtimes,1)    
    