function [] = on_sust1_Analysis(sigsToMeasure, rawSigsToMeasure,  tvals, cellIDs, approxStart, approxEnd, viewPlots, saveFlag, saveDir)

[~, approxStartInd] = min(abs(tvals-approxStart));
[~, approxEndInd] = min(abs(tvals-approxEnd));

tRiseOS1 = zeros(1,size(sigsToMeasure,1));
sigStart = zeros(1,size(sigsToMeasure,1));
sigEnd = zeros(1,size(sigsToMeasure,1));
riseStrtInt = zeros(1,size(sigsToMeasure,1));
riseEndInt = zeros(1,size(sigsToMeasure,1));
sigAmp = zeros(1,size(sigsToMeasure,1));

ampAvgWidth = 5;

figure;
for pp = 1:size(sigsToMeasure,1)
    tmpData = sigsToMeasure(pp,:);
    
    dfData = diff(tmpData); % take first derivative
    df2Data = diff(tmpData,2); % take first derivative
    
    % Find where first derivative crosses zero
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <=0);
    zx = zci(dfData);
    zx2 = zci(df2Data);
    
    dfTimes = linspace(tvals(2),tvals(end),length(dfData)); % Determine time points for first derivative
    df2Times = linspace(tvals(3),tvals(end),length(df2Data)); % Determine time points for first derivative
    
    % Identify peak of first derivative (signal onset)
    [~, dfPkInds] = findpeaks(dfData./max(dfData(:)),'MinPeakHeight',0.2,'MinPeakProminence',0.2); 
    dfPkLoc = (approxStart - dfTimes(dfPkInds));
    dfPkInds(dfPkLoc>0) = [];
    dfPkLoc(dfPkLoc>0) = [];
    dfPkLoc = abs(dfPkLoc);
    [~,dfPkTm] = min(dfPkLoc);
    dfPk = dfTimes(dfPkInds(dfPkTm));
    
    if ~isempty(dfPk)
        % locate troughs in first derivative 
        troughData = islocalmin(dfData);
        troughInds = 1:length(troughData);
        troughInds = troughInds(troughData);
        zx = [zx; troughInds'];
        zx = unique(zx);

        sigStartG = (dfPk - dfTimes(zx))';
        sigStartG(sigStartG<0) = inf;
        [~,startInd] = min(sigStartG);

        % Calculate beginning and end of rise
        sigEnd(pp) = dfTimes(zx(startInd+1));
        sigStart(pp) = tvals(approxStartInd);
        if sigStart(pp)<approxStart
            sigStart(pp) = approxStart;
        end
        tRiseOS1(pp) = (sigEnd(pp)-sigStart(pp))/2;
        
        sigAmp(pp) = max(tmpData(:)); %mean(tmpData(zx(startInd+1)-round(ampAvgWidth/2):zx(startInd+1)+round(ampAvgWidth/2)));
        
%         noiseSTD = std(abs(rawSignal-max(rawSignal(:)).*tmpData));
        
        % Calculate Rise Half-Time
        [halfTime, riseStrtInt(pp), riseEndInt(pp)] = calcHalfRise(tmpData(approxStartInd:zx(startInd+1)+1),tvals(approxStartInd:zx(startInd+1)+1));
        tRiseOS1(pp) = halfTime - approxStart;

        % Display Signal
        if viewPlots
            plot(tvals,tmpData,'LineWidth',2); % plot signal
            hold on;
%             plot(tvals,rawSignal./max(abs(rawSignal(:))));
%             plot(dfTimes,dfData./max(abs(dfData(:))),'LineWidth',2); % plot first derivative
%             plot(df2Times,df2Data./max(abs(df2Data(:))),'LineWidth',2);
            plot([sigStart(pp) sigStart(pp)], [-1.5 1.5],'--k','LineWidth',2);
            plot([tvals(1) tvals(end)],[sigAmp(pp) sigAmp(pp)],'--k','LineWidth',2);
%             plot([dfPk dfPk],[-1.5 1.5], '--r','LineWidth',2);
            plot([sigEnd(pp) sigEnd(pp)],[-1.5 1.5], '--k','LineWidth',2);
            plot([halfTime halfTime],[-1.5 1.5], '--g','LineWidth',2);
            ylim([-1.5 1.5]);
%             legend('Filtered Signal','First Derivative','Second Derivative'); legend box off;
            xlabel('Time (sec)');
            ylabel('Normalized Intensity');
            set(gca,'FontSize',14);
%             set(gcf,'Position',[-1153,196,926,652]);
            drawnow;
            hold off; 
        end
        
%         totalNoise(pp) = sum(abs(rawSignal./max(abs(rawSignal(:))) - tmpData));
    end
end

%% Plot Histograms

tsRiseOS1Cleaned = tRiseOS1;
sigEndCleaned = sigEnd;
riseIntDiff = riseEndInt-riseStrtInt;
OS1CleanedRaw = rawSigsToMeasure;
ampsCleaned = sigAmp;
cellIDsCleaned = cellIDs;

tsRiseOS1Cleaned(riseIntDiff<=0) = [];
sigEndCleaned(riseIntDiff<=0) =[];
cellIDsCleaned(riseIntDiff<=0) = [];
ampsCleaned(riseIntDiff<=0) = [];
OS1CleanedRaw(riseIntDiff<=0,:) = [];

tsRiseOS1Cleaned(tsRiseOS1Cleaned<=0) = [];
sigEndCleaned(tsRiseOS1Cleaned<=0) = [];
cellIDsCleaned(tsRiseOS1Cleaned<=0) = [];
ampsCleaned(tsRiseOS1Cleaned<=0) = [];
OS1CleanedRaw(tsRiseOS1Cleaned<=0,:) = [];

cellIDsCleaned = ['Time' cellIDsCleaned]; 

sigStats = [mean(tsRiseOS1Cleaned) std(tsRiseOS1Cleaned); mean(sigEndCleaned-approxStart) std(sigEndCleaned-approxStart); ...
    mean(ampsCleaned) std(ampsCleaned)];

% Display Histograms of Rise Half-Time
figure; histogram(tsRiseOS1Cleaned); 
    xlabel('Time (sec)'); 
    ylabel('Frequency'); 
    title('\tau_{1/2} Rise Times - ON Sustained 1');
    set(gca,'fontweight','bold','FontSize',14);
    if saveFlag
        saveas(gcf, [saveDir '/responseTime.tiff']);
    end
    
figure; histogram(ampsCleaned); 
    xlabel('Amplitude'); 
    ylabel('Frequency'); 
    title('Signal Amplitudes - ON Sustained 1');
    set(gca,'fontweight','bold','FontSize',14);
    if saveFlag
        saveas(gcf, [saveDir '/amplitudes.tiff']);
    end
    
%% Save quick look text file

fid = fopen([saveDir '\results.txt'],'w');
fprintf(fid,'Cell Type: ON-Sustained 1\n');
fprintf(fid,'N = %d\n', length(tsRiseOS1Cleaned));
fprintf(fid,'t_rise = %0.2f +/- %0.2f seconds\n',sigStats(1,1),sigStats(1,2));
fprintf(fid,'sigAmp = %0.2f +/- %0.2f\n',sigStats(3,1),sigStats(3,2));
fclose(fid);

%% Save Inidividual Cell Data

T = num2cell([tvals'; OS1CleanedRaw]');
T = table([cellIDsCleaned; T]);
tmpFname = 'rawAnalzyedSignals.csv';
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);

cellIDsCleaned{1} = 'Stats';
T = num2cell([tsRiseOS1Cleaned; ampsCleaned]);
T = table([cellIDsCleaned; [{'Rise Time (sec)' 'Signal Amplitude'}' T]]);
tmpFname = 'analyzedSignalStats.csv';
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);