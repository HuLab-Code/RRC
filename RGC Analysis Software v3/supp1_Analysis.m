function [] = supp2_Analysis(sigsToMeasure, rawSigsToMeasure, tvals, cellIDs, approxStart, approxEnd, viewPlots, saveFlag, saveDir)

[~, approxStartInd] = min(abs(tvals-approxStart));
[~, approxEndInd] = min(abs(tvals-approxEnd));

tsRiseOS2_1 = zeros(1,size(sigsToMeasure,1));
tsFallOS2_2 = zeros(1,size(sigsToMeasure,1));
sigStart1 = zeros(1,size(sigsToMeasure,1));
sigEnd1 = zeros(1,size(sigsToMeasure,1));
sigStart2 = zeros(1,size(sigsToMeasure,1));
sigEnd2 = zeros(1,size(sigsToMeasure,1));
riseStrtInt = zeros(1,size(sigsToMeasure,1));
riseEndInt = zeros(1,size(sigsToMeasure,1));
fallStrtInt = zeros(1,size(sigsToMeasure,1));
fallEndInt = zeros(1,size(sigsToMeasure,1));
sigWidth = zeros(1,size(sigsToMeasure,1));
sigAmp = zeros(1,size(sigsToMeasure,1));

ampAvgWidth = 5;

figure(1); 
for pp = 1:size(sigsToMeasure,1)
    tmpData = sigsToMeasure(pp,:);
    tmpData = -tmpData;
    
    dfData = diff(tmpData); % take first derivative
    df2Data = diff(tmpData,2);
    
    % Find where first derivative crosses zero
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <=0);
    zx = zci(dfData);
    zx2 = zci(df2Data);
    
    dfTimes = linspace(tvals(2),tvals(end),length(dfData)); % Determine time points for first derivative
    df2Times = linspace(tvals(3),tvals(end),length(df2Data)); 
    
    % Identify peak of first derivative (signal onset)
    [~, dfPkInds] = findpeaks(dfData./max(dfData(:)),'MinPeakHeight',0.2,'MinPeakProminence',0.2); 
    dfPkLoc = (approxStart - dfTimes(dfPkInds));
    dfPkInds(dfPkLoc>0) = [];
    dfPkLoc(dfPkLoc>0) = [];
    dfPkLoc = abs(dfPkLoc);
    [~,dfPkTm] = min(dfPkLoc);
    dfPk = dfTimes(dfPkInds(dfPkTm));
    
    [~, dfPkInds2] = findpeaks(-dfData./max(dfData(:)),'MinPeakHeight',0.2,'MinPeakProminence',0.2);  
    dfPkLoc2 = (approxEnd - dfTimes(dfPkInds2));
    dfPkInds2(dfPkLoc2>0) = [];
    dfPkLoc2(dfPkLoc2>0) = [];
    dfPkLoc2 = abs(dfPkLoc2);
    [~,dfPkTm2] = min(dfPkLoc2);
    dfPk2 = dfTimes(dfPkInds2(dfPkTm2));
    
    if ~isempty(dfPk) && ~isempty(dfPk2)
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
        sigEnd1(pp) = dfTimes(zx(startInd+1));
        sigStart1(pp) = approxStart;
        
        % Calculate Rise Half-Time
        [halfRiseTime, riseStrtInt(pp), riseEndInt(pp)] = calcHalfRise(tmpData(approxStartInd:approxEndInd),tvals(approxStartInd:approxEndInd));
        tsRiseOS2_1(pp) = halfRiseTime - approxStart;
        
        peakData = islocalmin(-dfData);
        peakInds = 1:length(peakData);
        peakInds = peakInds(peakData);
        zx = [zx; peakInds'];
        zx = unique(zx);
        
        sigStartG = (dfPk2 - dfTimes(zx))';
        sigStartG(sigStartG<=0) = inf;
        [~,startInd] = min(sigStartG);

        % Calculate beginning and end of rise
        sigEnd2(pp) = 50;
        sigStart2(pp) = approxEnd; %dfTimes(zx(startInd));
        sigAmp(pp) = mean(tmpData(zx(startInd)-round(ampAvgWidth/2):zx(startInd)+round(ampAvgWidth/2)));
        
        % Calculate Fall Half-Time
        [halfFallTime, fallStrtInt(pp), fallEndInt(pp)] = calcHalfDecay(tmpData(approxEndInd:236),tvals(approxEndInd:236));
        tsFallOS2_2(pp) = halfFallTime - sigStart2(pp);
        
        sigWidth(pp) = halfFallTime - halfRiseTime;

        % Display Signal
        if viewPlots
            plot(tvals,-tmpData,'LineWidth',2); % plot signal
            hold on;
%             plot(dfTimes,-dfData./max(abs(dfData(:))),'LineWidth',2); % plot first derivative
%             plot(df2Times,-df2Data./max(abs(df2Data(:))),'LineWidth',2); 
            plot([sigStart1(pp) sigStart1(pp)], [-1.5 1.5],'--k','LineWidth',2);
            plot([sigStart2(pp) sigStart2(pp)], [-1.5 1.5],'--k','LineWidth',2);
            plot([tvals(1) tvals(end)],[-sigAmp(pp) -sigAmp(pp)],'--k','LineWidth',2);
%             plot([dfPk dfPk],[-1.5 1.5], '--r','LineWidth',2);
%             plot([dfPk2 dfPk2],[-1.5 1.5], '--r','LineWidth',2);
            plot([halfRiseTime halfRiseTime],[-1.5 1.5], '--g','LineWidth',2);
            plot([halfFallTime halfFallTime],[-1.5 1.5], '--g','LineWidth',2);
%             plot([sigEnd1(pp) sigEnd1(pp)],[-1.5 1.5], '--k','LineWidth',2);
            plot([sigEnd2(pp) sigEnd2(pp)],[-1.5 1.5], '--k','LineWidth',2);
            ylim([-1.5 1.5]);
%             legend('Filtered Signal','First Derivative','Second Derivative'); legend box off;
            xlabel('Time (sec)');
            ylabel('Normalized Intensity');
            set(gca,'FontSize',14);
            set(gcf,'Position',[-1153,196,926,652]);
            drawnow;
            hold off; 
        end
    end
end

%%
riseIntDiff = riseEndInt-riseStrtInt;
fallIntDiff = fallStrtInt - fallEndInt;
tsRiseCleaned = tsRiseOS2_1;
tsFallCleaned = tsFallOS2_2;
sigWidthCleaned = sigWidth;
ampsCleaned = sigAmp;
offSuppCleaned = sigsToMeasure;
offSuppCleanedRaw = rawSigsToMeasure;
cellIDsCleaned = cellIDs;


% Remove points that don't have rising signal
tsRiseCleaned(riseIntDiff<=0) = [];
tsFallCleaned(riseIntDiff<=0) = [];
sigWidthCleaned(riseIntDiff<=0) = [];
ampsCleaned(riseIntDiff<=0) = [];
cellIDsCleaned(riseIntDiff<=0) = [];
offSuppCleaned(riseIntDiff<=0,:) = [];
offSuppCleanedRaw(riseIntDiff<=0,:) = [];


% Remove points that don't have falling signal
fallIntDiff(riseIntDiff<=0) = [];
tsRiseCleaned(fallIntDiff<=0) = [];
tsFallCleaned(fallIntDiff<=0) = [];
sigWidthCleaned(fallIntDiff<=0) = [];
ampsCleaned(fallIntDiff<=0) = [];
cellIDsCleaned(fallIntDiff<=0) = [];
offSuppCleaned(fallIntDiff<=0,:) = [];
offSuppCleanedRaw(fallIntDiff<=0,:) = [];

% Remove points that weren't analyzed
tsRiseCleaned(tsRiseCleaned<=0) = [];
sigWidthCleaned(tsFallCleaned<=0) = [];
ampsCleaned(tsFallCleaned<=0) = [];
cellIDsCleaned(tsFallCleaned<=0) = [];
offSuppCleaned(:,tsFallCleaned<=0) = [];
offSuppCleanedRaw(:,tsFallCleaned<=0) = [];
tsFallCleaned(tsFallCleaned<=0) = [];

cellIDsCleaned = ['Time' cellIDsCleaned]; 

sigStats = [mean(tsRiseCleaned) std(tsRiseCleaned); mean(tsFallCleaned) std(tsFallCleaned); mean(sigWidthCleaned) std(sigWidthCleaned);  mean(ampsCleaned) std(ampsCleaned)];

% Display Histograms of Rise/Fall Half-Time
figure; histogram(tsRiseCleaned,'binwidth',1); hold on; histogram(tsFallCleaned,'binwidth',1); 
    xlabel('Time (sec)');
    ylabel('Frequency');
    title('\tau_{1/2} - Suppression 2');
    legend('Rising','Falling'); legend boxoff;
    set(gca,'FontSize',14,'fontweight','bold');
    
% Display Signal Width Histogram
figure; histogram(sigWidthCleaned,'binwidth',1);
    xlabel('Time (sec)');
    ylabel('Frequency');
    title('Signal Width - Suppresion 2');
    set(gca,'FontSize',14,'fontweight','bold');
    
% Display Signal Width Histogram
figure; histogram(ampsCleaned);
    xlabel('Amplitude');
    ylabel('Frequency');
    title('Signal Amplitude - Suppression 2');
    set(gca,'FontSize',14,'fontweight','bold');
%%
fid = fopen([saveDir '\results.txt'],'w');
fprintf(fid,'Cell Type: Suppression 2\n');
fprintf(fid,'N = %d\n', length(tsFallCleaned));
fprintf(fid,'t_rise = %0.2f +/- %0.2f seconds\n',sigStats(1,1),sigStats(1,2));
fprintf(fid,'t_fall = %0.2f +/- %0.2f seconds\n',sigStats(2,1),sigStats(2,2));
fprintf(fid,'sigWidth = %0.2f +/- %0.2f seconds\n',sigStats(3,1),sigStats(3,2));
fprintf(fid,'sigAmp = %0.2f +/- %0.2f\n',sigStats(4,1),sigStats(4,2));
fclose(fid);

%%
T = num2cell([tvals'; offSuppCleanedRaw]');
T = table([cellIDsCleaned; T]);
tmpFname = 'rawAnalzyedSignals.csv';
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);

cellIDsCleaned{1} = 'Stats';
T = num2cell([tsRiseCleaned; tsFallCleaned; sigWidthCleaned; ampsCleaned]);
T = table([cellIDsCleaned; [{'Rise Time (sec)' 'Fall Time (sec)' 'Signal Width (sec)' 'Signal Amplitude'}' T]]);
tmpFname = 'analyzedSignalStats.csv';
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);
