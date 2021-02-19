%Pulse log processing script, step 1
%Author: Ali Gholamrezaei, University of Sydney

close all;
clear all;

%folder = ''; Specify path to .puls files
files = dir(strcat(folder, '*.puls'));

filter = 3;

if filter == 1
Fs = 1000; 
Order = 1000;       
Fc = [1 30];       
flag = 'scale';  
win = hamming(Order+1);
b  = fir1(Order, Fc/(Fs/2), 'band', win, flag);
Hd = dfilt.dffir(b);
D = mean(grpdelay(Hd)); 
elseif  filter == 2
Fs = 1000; 
Order = 1;  %first order
Fc = 1;  %should be 3dB down at the cutoff
[but,a] = butter(Order,Fc/(Fs/2),'high'); 
elseif  filter == 3
Fs = 1000; 
Order = 1;  %first order
Fc = [1 30];  %should be 3dB down at the cutoff
[but,a] = butter(Order,Fc/(Fs/2),'bandpass'); 
end

MinpulsePD = 700;
MinpulsePH = 50; %max(pulseraw)/div;

%errors:
minVP = 50;
maxVP = 200;
maxRIntdiff = 100;
maxIBIIntdiff = 100;
maxIBIs3 = 100;



sum_criteria_t = 5;

movW = 15;
movT = 350;

for filenr = 36%:size(files,1)

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

signal = alldata{1}(startsignal:endsignal);
signal = str2double(signal);

time = 1:1000/PULS_SAMPLES_PER_SECOND:length(signal)*1000/PULS_SAMPLES_PER_SECOND;
time = time';

% plot(time, signal);
% hold on

system_triggers = find(signal(:,1)>=5000);
Rtimes_system = 1000/PULS_SAMPLES_PER_SECOND * system_triggers;
IBIs_system = [NaN;diff(Rtimes_system)];

signal_edited = signal;

signal_edited(system_triggers) = NaN;

signal_edited = fillmissing(signal_edited, 'spline');  % or pchip

% signal_edited = [NaN; diff(signal_edited)];
% plot(time, signal_edited); 

time2 = 1:1000/Fs:time(end);
signal_int = interp1(time,signal_edited,time2,'spline')'; %


% yyaxis left
% plot(time2, signal_int);


if filter == 1
signal_filtered = filter(Hd, [signal_int; zeros(D,1)]); 
signal_filtered = signal_filtered(D+1:end,1); 
else 
signal_filtered = filtfilt(but,a,signal_int);
end

% yyaxis right
% plot(time2, signal_filtered);


% to save files with all columns; time original record, signal without system generated peaks, peaks, ...

pulseraw = signal_filtered; %signal_filtered;
% div = 2;
% Fs = fs; 
% N  = 1000; 
% Fc = 1; 
% flag = 'scale';  
% win = hamming(N+1);
% b  = fir1(N, Fc/(Fs/2), 'high', win, flag);
% Hd = dfilt.dffir(b);
% D = mean(grpdelay(Hd)); 

%             ECGraw = filter(Hd, [ECGraw(:,1); zeros(D,1)]); 
%             ECGraw = ECGraw(D+1:end); 

            [R_pks, Rtimes] = findpeaks(pulseraw, 'MinPeakHeight', MinpulsePH, 'MinPeakDistance', MinpulsePD);
            
            if Rtimes(1)<200
                Rtimes=Rtimes(2:end);
                R_pks=R_pks(2:end);
            end
            
            IBIs = [NaN; diff(Rtimes)];
            
            SD = [NaN; diff(IBIs)];
            
            SDab = abs(SD);
            
            SDm = movmean(SDab, movW);
            SDmdiff = 100*abs(SDm-SDab)./SDm;

            V = [NaN];
            I = [NaN];
            VP = [NaN];

            
            
            for   beatNr = 1:size(Rtimes,1)-1
                v = [];
                i = [];
                [v, i] = min(pulseraw(Rtimes(beatNr):Rtimes(beatNr+1)));
                i = Rtimes(beatNr) + i;
                vp = Rtimes(beatNr+1) - i;
                
                V = [V; v];
                I = [I; i];
                VP = [VP; vp];
            end

            RIntdiff = [];
            IBIIntdiff = [];
            IBIs3 = [];
            Rtimes2 = Rtimes;
            IBIs2 = IBIs;
            for beatNr = 1: size(Rtimes,1)
            
%                 Rtimes2 = Rtimes;
                Rtimes2(beatNr)=NaN;
                Rtimes2 = fillmissing(Rtimes2, 'spline');  % or pchip
                RIntdiff(beatNr) = abs(Rtimes(beatNr) - Rtimes2(beatNr));
                
                if beatNr > 2 & beatNr < size(Rtimes,1)
                IBIs2(beatNr)=NaN;
                IBIs2 = fillmissing(IBIs2, 'spline');  % or pchip
                IBIIntdiff(beatNr) = abs(IBIs(beatNr) - IBIs2(beatNr));
                
                IBIs3(beatNr) = abs(IBIs(beatNr) - mean(IBIs(beatNr-1),IBIs(beatNr+1)));
                else
                    IBIs3(beatNr) = NaN;
                end                
            
            end

            IBIs2(1) = NaN;
            IBIs3=IBIs3';
%             IBIs3 = [NaN;NaN;IBIs3];
            RIntdiff = RIntdiff';
            IBIIntdiff = IBIIntdiff';
            IBIIntdiff = [NaN;IBIIntdiff];
            
            sum_criteria =[];
            
%             RIntdiffsort = sort(RIntdiff, 'descend');
%             percentile5 = prctile(RIntdiff,[1:100], 'all'); 
            
%             L = log(RIntdiff);
            
            errors = []; 

            for   beatNr = 1:size(Rtimes,1)-1
              
            if (RIntdiff(beatNr)>maxRIntdiff & VP(beatNr)>maxVP) ...
                   | VP(beatNr)<minVP | (VP(beatNr)>maxVP & IBIIntdiff>maxIBIIntdiff & IBIs3>maxIBIs3) ...
                   | I(beatNr+1)-Rtimes(beatNr)<100 % | V(beatNr)<-877 % | SDmdiff(beatNr+1)>movT 
                
                errors(beatNr)=1;
            else 
                errors(beatNr)=0;
                
            end
            
            
%             sum_criteria(beatNr) = RIntdiff(beatNr)>maxRIntdiff + VP(beatNr)>maxVP + ... 
%                 VP(beatNr)<minVP + IBIs3(beatNr)>maxIBIs3 + SDmdiff(beatNr)>movT;
            
            end

            errors = [errors'; NaN];
            
%             sum_criteria =sum_criteria';
%             sum_criteria = [NaN;sum_criteria];
            
            SDmdiff2 = [SDmdiff(2:end);NaN];

            sum_criteria = ((RIntdiff>maxRIntdiff)*1) + (VP>maxVP) + (VP>400) + (VP>700) + ...
                ((IBIIntdiff>maxIBIIntdiff)*1) + (IBIIntdiff>250) + ... 
                (VP<minVP) + ((IBIs3>maxIBIs3)*1) + (IBIs3>300) + (IBIs3>500) + ...
                (SDmdiff2>movT);
            
            Time = 1:size(pulseraw,1);
                Time = downsample(Time, 5);
                pulseraw = downsample(pulseraw, 5);
              
                nf = 1;
%%
            figure(nf);
            ax1 = subplot(4,1,1);
            yyaxis left
            plot(Time, pulseraw);
            hold on
            plot(Rtimes, R_pks, 'o', 'Color', 'b');
            hold on
            plot(I, V, 'o', 'Color', 'r');
            
            if nanmean(errors)~=0
            hold on
            line([Rtimes(errors==1,1) Rtimes(errors==1,1)]', ylim, 'Marker','none','LineStyle','-','Color',[1 0 0]);
            end
            if nanmean(sum_criteria)~=0 & isempty(Rtimes(sum_criteria>=sum_criteria_t,1))~=1
%             hold on
            line([Rtimes(sum_criteria>=sum_criteria_t,1) Rtimes(sum_criteria>=sum_criteria_t,1)]', ylim, ... 
                'Marker','none','LineStyle','--','Color', 'K');
            end
            ylim([-2000 2000]);
            xlim([min(Time), max(Time)]);
            title('pulse');

            yyaxis right
            plot(time, signal, 'LineStyle', '--');
            
            ax2 = subplot(4,1,2);
            yyaxis left 
            plot(Rtimes, IBIs, '-o');
            hold on
            plot(Rtimes, IBIs2, '-*');
            hold on
            plot(Rtimes_system, IBIs_system, '-+', 'Color', 'K');
            
            
%             ref1500 = refline(0 , 1500);
%             ref500 = refline(0 , 500);
            ref3std = refline(0 , (nanmean(IBIs) + 3*nanstd(IBIs)));
            ref3std.Color = 'b';
            ref3_std = refline(0 , (nanmean(IBIs) - 3*nanstd(IBIs)));
            ref3_std.Color = 'b';
            if nanmean(errors)~=0
            hold on
            line([Rtimes(errors==1,1) Rtimes(errors==1,1)]', ylim, 'Marker','none','LineStyle','-','Color',[1 0 0]);
            end
            if nanmean(sum_criteria)~=0 & isempty(Rtimes(sum_criteria>=sum_criteria_t,1))~=1
%             hold on
            line([Rtimes(sum_criteria>=sum_criteria_t,1) Rtimes(sum_criteria>=sum_criteria_t,1)]', ylim, ... 
                'Marker','none','LineStyle','--','Color', 'K');
            end                       
            ylim([min(IBIs)-100, max(IBIs)+100]);
            xlim([min(Time), max(Time)]);
%             title('IBI');
            ylabel('msec');
                        
            yyaxis right
            plot(Rtimes, SD, '*');
            hold on
            plot(Rtimes, SDmdiff, '-o');
            refline(0 , movT);
            ylim([min(min([SD, SDmdiff]))-100, max(max([SD, SDmdiff]))+100]);
            xlim([min(Time), max(Time)]);
            title('IBI, SD, SD diff from movemean');
            ylabel('msec / % from movmean');
     
            ax3 = subplot(4,1,3);
            yyaxis left
            plot(Rtimes, VP, 'o');
            refline(0 , 200);
            ylabel('msec');
            if nanmean(errors)~=0
            hold on
            line([Rtimes(errors==1,1) Rtimes(errors==1,1)]', ylim, 'Marker','none','LineStyle','-','Color',[1 0 0]);
            
            
            end
            if nanmean(sum_criteria)~=0 & isempty(Rtimes(sum_criteria>=sum_criteria_t,1))~=1
%             hold on
            line([Rtimes(sum_criteria>=sum_criteria_t,1) Rtimes(sum_criteria>=sum_criteria_t,1)]', ylim, ... 
                'Marker','none','LineStyle','--','Color', 'K');
            
            
            end
            yyaxis right
            plot(Rtimes, RIntdiff, '*');
            
            hold on
            plot(Rtimes, IBIIntdiff, '^');
            hold on
            plot(Rtimes, IBIs3, '+', 'Color', 'K');
            
            
            refline(0 , 100);
%             refline(0, (nanmean(RIntdiff)+(3*nanstd(RIntdiff))));
            title('D2S time, IBI diff if interpolated');
            ylabel('msec');
            
            ax4 = subplot(4,1,4);
            yyaxis left
            plot(Rtimes, sum_criteria, '-o');

%             refline(0 , 200);
%             ylabel('msec');

            linkaxes([ax1,ax2,ax3, ax4], 'x');
%%

end

    errors2 = errors;
    errors2(sum_criteria>=sum_criteria_t,1)=1;
    errors2(sum_criteria<sum_criteria_t,1)=0;

	data = [Rtimes R_pks I V VP IBIs RIntdiff IBIIntdiff ...
            IBIs3 SDmdiff2 sum_criteria errors2 errors];
        
        data = round(data, 1);
        
    data = array2table(data, 'variablenames', {'Rtimes' 'Pks' 'Vtimes' 'Vals' 'VPtimes' ...
        'IBIs' 'Rtimes_Int_diff' 'IBIs_Int_diff' 'IBIs_adj_diff' ...
        'SD_mov_diff' 'sum_criteria' 'sum_criteria_th' 'errors'});
    
    writetable(data, strcat(fname, '.txt'), 'Delimiter', '\t');
    
    saveas(nf, strcat(fname, '_.fig'));
    
    rows = find(sum_criteria>=sum_criteria_t);
    
    data(rows, :)
    
    (sum(errors2==1))*100/size(IBIs,1)
    (sum(errors2==1)+sum(errors==1))*100/size(IBIs,1)
% fclose(fname);