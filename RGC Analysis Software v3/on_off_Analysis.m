function [] = on_off_Analysis(sigsToMeasure, rawSigsToMeasure,  tvals, cellIDs, approxStart, approxEnd, viewPlots, saveFlag, saveDir)

[~, approxStartInd] = min(abs(tvals-approxStart));
[~, approxEndInd] = min(abs(tvals-approxEnd));

tsRiseOnOff_1 = zeros(1,size(sigsToMeasure,1));
tsRiseOnOff_2 = zeros(1,size(sigsToMeasure,1));
tsFallOnOff_1 = zeros(1,size(sigsToMeasure,1));
tsFallOnOff_2 = zeros(1,size(sigsToMeasure,1));
sigWidth1 = zeros(1,size(sigsToMeasure,1));
sigWidth2 = zeros(1,size(sigsToMeasure,1));
sigStart1 = zeros(1,size(sigsToMeasure,1));
sigEnd1 = zeros(1,size(sigsToMeasure,1));
sigStart2 = zeros(1,size(sigsToMeasure,1));
sigEnd2 = zeros(1,size(sigsToMeasure,1));
sigEnd3 = zeros(1,size(sigsToMeasure,1));
rise1StrtInt = zeros(1,size(sigsToMeasure,1));
rise2StrtInt = zeros(1,size(sigsToMeasure,1));
rise1EndInt = zeros(1,size(sigsToMeasure,1));
rise2EndInt = zeros(1,size(sigsToMeasure,1));
fall1StrtInt = zeros(1,size(sigsToMeasure,1));
fall2StrtInt = zeros(1,size(sigsToMeasure,1));
fall1EndInt = zeros(1,size(sigsToMeasure,1));
fall2EndInt = zeros(1,size(sigsToMeasure,1));
sigAmp1 = zeros(1,size(sigsToMeasure,1));
sigAmp2 = zeros(1,size(sigsToMeasure,1));
sigAmpMean = zeros(1,size(sigsToMeasure,1));

ampAvgWidth = 5;

figure; 
for pp = 1:size(sigsToMeasure,1)
    tmpData = sigsToMeasure(pp,:);
    
    dfData = diff(tmpData); % take first derivative
    df2Data = diff(tmpData,2);
    
    % Find where first derivative crosses zero
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <=0);
    zx = zci(dfData);
    z2x = zci(df2Data);
    
    dfTimes = linspace(tvals(2),tvals(end),length(dfData)); % Determine time points for first derivative
    df2Times = linspace(tvals(3),tvals(end),length(df2Data));
    
    % Identify peak of first derivative (signal onset)
    [~, dfPkInds] = findpeaks(dfData./max(dfData(:)),'MinPeakHeight',0.15,'MinPeakProminence',0.2); 
    dfPkLoc = (approxStart - dfTimes(dfPkInds));
    dfPkInds(dfPkLoc>0) = [];
    dfPkLoc(dfPkLoc>0) = [];
    dfPkLoc = abs(dfPkLoc);
    [~,dfPkTm] = min(dfPkLoc);
    dfPk = dfTimes(dfPkInds(dfPkTm));
    if dfPk > approxStart + 10
        dfPk = [];
    end
    
    [~, dfPkInds2] = findpeaks(dfData./max(dfData(:)),'MinPeakHeight',0.15,'MinPeakProminence',0.2);  
    dfPkLoc2 = (approxEnd - dfTimes(dfPkInds2));
    dfPkInds2(dfPkLoc2>0) = [];
    dfPkLoc2(dfPkLoc2>0) = [];
    dfPkLoc2 = abs(dfPkLoc2);
    [~,dfPkTm2] = min(dfPkLoc2);
    dfPk2 = dfTimes(dfPkInds2(dfPkTm2));
    
    sigTrough = islocalmin(tmpData);
    sigTroughInds = 1:length(sigTrough);
    sigTroughInds = sigTroughInds(sigTrough);
    
    if ~isempty(dfPk) && ~isempty(dfPk2)
%         filtSigsMeasured(:,pp) = tmpData;
        % locate troughs in first derivative 
        troughData = islocalmin(dfData);
        troughInds = 1:length(troughData);
        troughInds = troughInds(troughData);
        zxf = [zx; troughInds'];
        zxf = unique(zxf);

        sigStartG = (dfPk - dfTimes(zxf))';
        sigStartG(sigStartG<0) = inf;
        [~,startInd] = min(sigStartG);

        % Calculate beginning and end of rise
        sigEnd1(pp) = dfTimes(zxf(startInd+1));
        sigStart1(pp) = approxStart;
        sigAmp1(pp) = mean(tmpData(zxf(startInd+1)-round(ampAvgWidth/2):zxf(startInd+1)+round(ampAvgWidth/2)));
        
        [halfRiseTime1, rise1StrtInt(pp), rise1EndInt(pp)] = calcHalfRise(tmpData(approxStartInd:zxf(startInd+1)+1),tvals(approxStartInd:zxf(startInd+1)+1));
        tsRiseOnOff_1(pp) = halfRiseTime1 - approxStart;
        
        peakData = islocalmin(-dfData);
        peakInds = 1:length(peakData);
        peakInds = peakInds(peakData);
        zxf2 = [zxf; peakInds'];
        zxf2 = unique(zxf2);
        
        sigStartG = (dfPk2 - dfTimes(zxf2))';
        sigStartG(sigStartG<=0) = inf;
        [~,startInd2] = min(sigStartG);
        try
            [halfFallTime1, fall1StrtInt(pp), fall1EndInt(pp)] = calcHalfDecay(tmpData(zxf(startInd+1)+1:approxEndInd),tvals(zxf(startInd+1)+1:approxEndInd));
            tsFallOnOff_1(pp) = halfFallTime1 - sigEnd1(pp);
        end

        % Calculate beginning and end of second rise
        [~, endInd1] = findpeaks(tmpData(approxEndInd:end),'MinPeakProminence', 0.2, 'MinPeakHeight',0);
        if length(endInd1) > 1
            tmp = tmpData(approxEndInd:end);
            [~,ss] = max(tmp(endInd1));
            endInd1 = endInd1(ss);
        end
        endInd1 = endInd1+approxEndInd-1;
        sigEndG = (dfPk2 - dfTimes(zx))';
        sigEndG(sigEndG>=0) = inf;
        [~,sigEndInd] = min(abs(sigEndG));
        
        if tmpData(endInd1) > tmpData(zx(sigEndInd))
            sigEndInd = endInd1;
        else
            sigEndInd = zx(sigEndInd);
        end
        
        sigEnd2(pp) = dfTimes(sigEndInd);
        sigStart2(pp) = approxEnd;
        sigAmp2(pp) = mean(tmpData(sigEndInd-round(ampAvgWidth/2):sigEndInd+round(ampAvgWidth/2)));
        
        [halfRiseTime2, rise2StrtInt(pp), rise2EndInt(pp)] = calcHalfRise(tmpData(approxEndInd:sigEndInd+1),tvals(approxEndInd:sigEndInd+1));
        tsRiseOnOff_2(pp) = halfRiseTime2 - approxEnd;
        
        sigEndG = (sigEnd2(pp) - tvals(sigTroughInds));
        sigEndG(sigEndG >=0) = inf;
        [~,sigEndInd2] = min(abs(sigEndG));
        sigEnd3(pp) = tvals(sigTroughInds(sigEndInd2));
        try
            [halfFallTime2, fall2StrtInt(pp), fall2EndInt(pp)] = calcHalfDecay(tmpData(sigEndInd+1:sigTroughInds(sigEndInd2)),tvals(sigEndInd+1:sigTroughInds(sigEndInd2)));
            tsFallOnOff_2(pp) = halfFallTime2 - sigEnd2(pp);
            
            sigWidth1(pp) = halfFallTime1 - halfRiseTime1;
            sigWidth2(pp) = halfFallTime2 - halfRiseTime2;
        
        % Display Signal
        if viewPlots
            plot(tvals,tmpData,'LineWidth',2); % plot signal
            hold on;
%             plot(tvals, rawSignal./max(abs(rawSignal(:))));
%             plot(dfTimes,dfData./max(abs(dfData(:))),'LineWidth',2); % plot first derivative
%             plot(df2Times,df2Data./max(abs(df2Data(:))),'LineWidth',2); 
            plot([sigStart1(pp) sigStart1(pp)], [-1.5 1.5],'--k','LineWidth',2);
            plot([sigStart2(pp) sigStart2(pp)], [-1.5 1.5],'--k','LineWidth',2);
            plot([tvals(1) tvals(end)],[sigAmp1(pp) sigAmp1(pp)],'--k','LineWidth',2);
            plot([tvals(1) tvals(end)],[sigAmp2(pp) sigAmp2(pp)],'--k','LineWidth',2);
%             plot([dfPk dfPk],[-1.5 1.5], '--r','LineWidth',2);
            plot([halfRiseTime1 halfRiseTime1],[-1.5 1.5], '--g','LineWidth',2);
            plot([halfFallTime1 halfFallTime1],[-1.5 1.5], '--g','LineWidth',2);
            plot([halfRiseTime2 halfRiseTime2],[-1.5 1.5], '--g','LineWidth',2);
            plot([halfFallTime2 halfFallTime2],[-1.5 1.5], '--g','LineWidth',2);
            plot([sigEnd1(pp) sigEnd1(pp)],[-1.5 1.5], '--k','LineWidth',2);
            plot([sigEnd2(pp) sigEnd2(pp)],[-1.5 1.5], '--k','LineWidth',2);
            plot([sigEnd3(pp) sigEnd3(pp)],[-1.5 1.5], '--k','LineWidth',2);
            ylim([-1.5 1.5]);
%             legend('Filtered Signal','First Derivative','Second Derivative'); legend box off;
            xlabel('Time (sec)');
            ylabel('Normalized Intensity');
            set(gca,'FontSize',14);
            set(gcf,'Position',[205,138,926,652]);
            drawnow;
            hold off; 
        end
        end
%         totalNoise(pp) = sum(abs(rawSignal./max(rawSignal(:)) - tmpData));
    end
end

%% Display Histograms of Rise Half-Time

sigAmpMean = (sigAmp1 + sigAmp2)/2;
rise1IntDiff = rise1EndInt-rise1StrtInt;
rise2IntDiff = rise2EndInt-rise2StrtInt;
fall1IntDiff = fall1StrtInt - fall1EndInt;
fall2IntDiff = fall2StrtInt - fall2EndInt;
tsRise1Cleaned = tsRiseOnOff_1;
tsRise2Cleaned = tsRiseOnOff_2;
tsFall1Cleaned = tsFallOnOff_1;
tsFall2Cleaned = tsFallOnOff_2;
sigWidth1Cleaned = sigWidth1;
sigWidth2Cleaned = sigWidth2;
sigAmp1Clean = sigAmp1;
sigAmp2Clean = sigAmp2;
sigAmpMeanClean = sigAmpMean;
onOffCleaned = sigsToMeasure;
onOffCleanedRaw = rawSigsToMeasure;
cellIDsCleaned = cellIDs;

% Remove points that don't have rising signal
tsRise1Cleaned(rise1IntDiff<=0) = [];
tsRise2Cleaned(rise1IntDiff<=0) = [];
tsFall1Cleaned(rise1IntDiff<=0) = [];
tsFall2Cleaned(rise1IntDiff<=0) = [];
sigWidth1Cleaned(rise1IntDiff<=0) = [];
sigWidth2Cleaned(rise1IntDiff<=0) = [];
sigAmp1Clean(rise1IntDiff<=0) = [];
sigAmp2Clean(rise1IntDiff<=0) = [];
sigAmpMeanClean(rise1IntDiff<=0) = [];
cellIDsCleaned(rise1IntDiff<=0) = [];
onOffCleaned(rise1IntDiff<=0,:) = [];
onOffCleanedRaw(rise1IntDiff<=0,:) = [];

% Remove points that don't have rising signal
rise2IntDiff(rise1IntDiff<=0) = [];
fall1IntDiff(rise1IntDiff<=0) = [];
fall2IntDiff(rise1IntDiff<=0) = [];
tsRise1Cleaned(rise2IntDiff<=0) = [];
tsRise2Cleaned(rise2IntDiff<=0) = [];
tsFall1Cleaned(rise2IntDiff<=0) = [];
tsFall2Cleaned(rise2IntDiff<=0) = [];
sigWidth1Cleaned(rise2IntDiff<=0) = [];
sigWidth2Cleaned(rise2IntDiff<=0) = [];
sigAmp1Clean(rise2IntDiff<=0) = [];
sigAmp2Clean(rise2IntDiff<=0) = [];
sigAmpMeanClean(rise2IntDiff<=0) = [];
cellIDsCleaned(rise2IntDiff<=0) = [];
onOffCleaned(rise2IntDiff<=0,:) = [];
onOffCleanedRaw(rise2IntDiff<=0,:) = [];

% Remove points that don't have falling signal
fall1IntDiff(rise2IntDiff<=0) = [];
fall2IntDiff(rise2IntDiff<=0) = [];
tsRise1Cleaned(fall1IntDiff<=0) = [];
tsRise2Cleaned(fall1IntDiff<=0) = [];
tsFall1Cleaned(fall1IntDiff<=0) = [];
tsFall2Cleaned(fall1IntDiff<=0) = [];
sigWidth1Cleaned(fall1IntDiff<=0) = [];
sigWidth2Cleaned(fall1IntDiff<=0) = [];
sigAmp1Clean(fall1IntDiff<=0) = [];
sigAmp2Clean(fall1IntDiff<=0) = [];
sigAmpMeanClean(fall1IntDiff<=0) = [];
cellIDsCleaned(fall1IntDiff<=0) = [];
onOffCleaned(fall1IntDiff<=0,:) = [];
onOffCleanedRaw(fall1IntDiff<=0,:) = [];

% Remove points that don't have falling signal
fall2IntDiff(fall1IntDiff<=0) = [];
tsRise1Cleaned(fall2IntDiff<=0) = [];
tsRise2Cleaned(fall2IntDiff<=0) = [];
tsFall1Cleaned(fall2IntDiff<=0) = [];
tsFall2Cleaned(fall2IntDiff<=0) = [];
sigWidth1Cleaned(fall2IntDiff<=0) = [];
sigWidth2Cleaned(fall2IntDiff<=0) = [];
sigAmp1Clean(fall2IntDiff<=0) = [];
sigAmp2Clean(fall2IntDiff<=0) = [];
sigAmpMeanClean(fall2IntDiff<=0) = [];
cellIDsCleaned(fall2IntDiff<=0) = [];
onOffCleaned(fall2IntDiff<=0,:) = [];
onOffCleanedRaw(fall2IntDiff<=0,:) = [];


% Remove points that weren't analyzed
tsRise1Cleaned(tsRise1Cleaned<=0) = [];
tsRise2Cleaned(tsRise2Cleaned<=0) = [];
sigWidth1Cleaned(tsFall2Cleaned<=0) = [];
sigWidth2Cleaned(tsRise2Cleaned<=0) = [];
sigAmp1Clean(tsRise2Cleaned<=0) = [];
sigAmp2Clean(tsRise2Cleaned<=0) = [];
sigAmpMeanClean(tsRise2Cleaned<=0) = [];
cellIDsCleaned(tsFall2Cleaned<=0) = [];
onOffCleaned(tsFall2Cleaned<=0,:) = [];
onOffCleanedRaw(tsFall2Cleaned<=0,:) = [];
tsFall1Cleaned(tsFall2Cleaned<=0) = [];
tsFall2Cleaned(tsFall2Cleaned<=0) = [];

cellIDsCleaned = ['Time' cellIDsCleaned]; 

sigStats = [mean(tsRise1Cleaned) std(tsRise1Cleaned); mean(tsFall1Cleaned) std(tsFall1Cleaned); mean(sigWidth1Cleaned) std(sigWidth1Cleaned); ...
    mean(tsRise2Cleaned) std(tsRise2Cleaned); mean(tsFall2Cleaned) std(tsFall2Cleaned); mean(sigWidth2Cleaned) std(sigWidth2Cleaned);
    mean(sigAmp1Clean) std(sigAmp1Clean); mean(sigAmp2Clean) std(sigAmp2Clean); mean(sigAmpMeanClean) std(sigAmpMeanClean)];

% Display Histograms of Rise/Fall Half-Time
figure; histogram(tsRise1Cleaned,'BinWidth',1); hold on; histogram(tsFall1Cleaned,'BinWidth',1); 
    xlabel('Time (sec)');
    ylabel('Frequency');
    title('\tau_{1/2} - ON-OFF Peak 1');
    legend('Rising','Falling'); legend boxoff;
    set(gca,'fontweight','bold','FontSize',14);
    if saveFlag
        saveas(gcf, [saveDir '/peak1ResponseTimes.tiff']);
    end
    
% Display Histograms of Rise/Fall Half-Time
figure; histogram(tsRise2Cleaned,'BinWidth',1); hold on; histogram(tsFall2Cleaned,'BinWidth',1); 
    xlabel('Time (sec)');
    ylabel('Frequency');
    title('\tau_{1/2} - ON-OFF Peak 2');
    legend('Rising','Falling'); legend boxoff;
    set(gca,'fontweight','bold','FontSize',14);
    if saveFlag
        saveas(gcf, [saveDir '/peak2ResponseTimes.tiff']);
    end
    
% Display Signal Width Histogram
figure; histogram(sigWidth1Cleaned,'BinWidth',1); hold on; histogram(sigWidth2Cleaned,'BinWidth',1);
    xlabel('Time (sec)');
    ylabel('Frequency');
    title('Signal Widths - ON-OFF');
    legend('Peak 1','Peak 2'); legend boxoff;
    set(gca,'fontweight','bold','FontSize',14);
    if saveFlag
        saveas(gcf, [saveDir '/peakWidths.tiff']);
    end
    
% Display Signal Width Histogram
figure; histogram(sigAmp1Clean); hold on; histogram(sigAmp2Clean);
    xlabel('Signal Amplitude');
    ylabel('Frequency');
    title('Signal Amplitude - ON-OFF');
    legend('Peak 1','Peak 2'); legend boxoff;
    set(gca,'fontweight','bold','FontSize',14);
    if saveFlag
        saveas(gcf, [saveDir '/peakAmps.tiff']);
    end
    
% Display Signal Width Histogram
figure; histogram(sigAmpMeanClean);
    xlabel('Signal Amplitude');
    ylabel('Frequency');
    title('Mean Signal Amplitude - ON-OFF');
    set(gca,'fontweight','bold','FontSize',14);
    if saveFlag
        saveas(gcf, [saveDir '/meanPeakAmps.tiff']);
    end

%% Save Quick Look Stats
fid = fopen([saveDir '\results.txt'],'w');
fprintf(fid,'Cell Type: ON-OFF\n');
fprintf(fid,'N = %d\n', length(tsFall2Cleaned));
fprintf(fid,'t_rise1 = %0.2f +/- %0.2f seconds\n',sigStats(1,1),sigStats(1,2));
fprintf(fid,'t_fall1 = %0.2f +/- %0.2f seconds\n',sigStats(2,1),sigStats(2,2));
fprintf(fid,'sigWidth1 = %0.2f +/- %0.2f seconds\n',sigStats(3,1),sigStats(3,2));
fprintf(fid,'t_rise2 = %0.2f +/- %0.2f seconds\n',sigStats(4,1),sigStats(4,2));
fprintf(fid,'t_fall2 = %0.2f +/- %0.2f seconds\n',sigStats(5,1),sigStats(5,2));
fprintf(fid,'sigWidth2 = %0.2f +/- %0.2f seconds\n',sigStats(6,1),sigStats(6,2));
fprintf(fid,'sigAmp1 = %0.2f +/- %0.2f \n',sigStats(7,1),sigStats(7,2));
fprintf(fid,'sigAmp2 = %0.2f +/- %0.2f \n',sigStats(8,1),sigStats(8,2));
fprintf(fid,'sigAmpMean = %0.2f +/- %0.2f \n',sigStats(9,1),sigStats(9,2));
fclose(fid);

%% Save Inidividual Cell Data

T = num2cell([tvals'; onOffCleanedRaw]');
T = table([cellIDsCleaned; T]);
tmpFname = 'rawAnalzyedSignals.csv'; 
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);

cellIDsCleaned{1} = 'Stats';
T = num2cell([tsRise1Cleaned; tsFall1Cleaned; sigWidth1Cleaned; tsRise2Cleaned; tsFall2Cleaned; sigWidth2Cleaned; sigAmp1Clean; sigAmp2Clean; sigAmpMeanClean]);
T = table([cellIDsCleaned; [{'Peak 1 Rise Time (sec)' 'Peak 1 Fall Time (sec)' 'Peak 1 Width (sec)' 'Peak 2 Rise Time (sec)' 'Peak 2 Fall Time (sec)' 'Peak 2 Width (sec)' 'Peak 1 Amplitude' 'Peak 2 Amplitude' 'Mean Amplitude'}' T]]);
tmpFname = 'analyzedSignalStats.csv';
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);