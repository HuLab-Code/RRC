function [] = on_trans_Analysis(sigsToMeasure, rawSigsToMeasure, tvals, cellIDs, approxStart, approxEnd, viewPlots, saveFlag, saveDir)

[~, approxStartInd] = min(abs(tvals-approxStart));
[~, approxEndInd] = min(abs(tvals-approxEnd));

tsRiseOT2_1 = zeros(1,size(sigsToMeasure,1));
tsFallOT2_2 = zeros(1,size(sigsToMeasure,1));
sigStart1 = zeros(1,size(sigsToMeasure,1));
riseStrtInt = zeros(1,size(sigsToMeasure,1));
riseEndInt = zeros(1,size(sigsToMeasure,1));
fallStrtInt = zeros(1,size(sigsToMeasure,1));
fallEndInt = zeros(1,size(sigsToMeasure,1));
sigEnd1 = zeros(1,size(sigsToMeasure,1));
sigStart2 = zeros(1,size(sigsToMeasure,1));
sigEnd2 = zeros(1,size(sigsToMeasure,1));
sigWidth = zeros(1,size(sigsToMeasure,1));
sigAmp = zeros(1,size(sigsToMeasure,1));

ampAvgWidth = 5;

figure(1); 
for pp = 1:size(sigsToMeasure,1)

    tmpData = sigsToMeasure(pp,:);
    
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
    
    sigTrough = islocalmin(tmpData);
    sigTroughInds = 1:length(sigTrough);
    sigTroughInds = sigTroughInds(sigTrough);
    
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
        sigEnd1(pp) = dfTimes(zx(startInd+1));
        sigStart1(pp) = approxStart;
        sigAmp(pp) = mean(tmpData(zx(startInd+1)-round(ampAvgWidth/2):zx(startInd+1)+round(ampAvgWidth/2)));
        
        if approxStartInd < (zx(startInd+1)+1)
            [halfRiseTime, riseStrtInt(pp), riseEndInt(pp)] = calcHalfRise(tmpData(approxStartInd:zx(startInd+1)+1),tvals(approxStartInd:zx(startInd+1)+1));
            tsRiseOT2_1(pp) = halfRiseTime - approxStart;

            peakData = islocalmin(-dfData);
            peakInds = 1:length(peakData);
            peakInds = peakInds(peakData);
    %         zx = [zx; peakInds'];
    %         zx = unique(zx);

            sigEndG = (approxEnd - tvals(sigTroughInds));
            sigEndG(sigEndG>0) = inf;
            [~,sigEndInd] = min(abs(sigEndG));

            % Calculate beginning and end of fall
            sigEnd2(pp) = tvals(sigTroughInds(sigEndInd));
            sigStart2(pp) = sigEnd1(pp);
            if (zx(startInd+1)+1) < sigTroughInds(sigEndInd)
                [halfFallTime, fallStrtInt(pp), fallEndInt(pp)] = calcHalfDecay(tmpData(zx(startInd+1)+1:sigTroughInds(sigEndInd)),tvals(zx(startInd+1)+1:sigTroughInds(sigEndInd)));
                tsFallOT2_2(pp) = halfFallTime - sigStart2(pp);
                
                sigWidth(pp) = halfFallTime - halfRiseTime;

                % Display Signal
                if viewPlots
                    plot(tvals,tmpData,'LineWidth',2); % plot signal
                    hold on;
%                     plot(tvals, rawSignal./max(abs(rawSignal(:))));
%                     plot(dfTimes,dfData./max(abs(dfData(:))),'LineWidth',2); % plot first derivative
%                     plot(df2Times,df2Data./max(abs(df2Data(:))),'LineWidth',2); 
                    plot([sigStart1(pp) sigStart1(pp)], [-1.5 1.5],'--k','LineWidth',2);
                    plot([sigStart2(pp) sigStart2(pp)], [-1.5 1.5],'--k','LineWidth',2);
                    plot([tvals(1) tvals(end)],[sigAmp(pp) sigAmp(pp)],'--k','LineWidth',2);
        %             plot([dfPk dfPk],[-1.5 1.5], '--r','LineWidth',2);
                    plot([halfRiseTime halfRiseTime],[-1.5 1.5], '--g','LineWidth',2);
                    plot([halfFallTime halfFallTime],[-1.5 1.5], '--g','LineWidth',2);
                    plot([sigEnd1(pp) sigEnd1(pp)],[-1.5 1.5], '--k','LineWidth',2);
                    plot([sigEnd2(pp) sigEnd2(pp)],[-1.5 1.5], '--k','LineWidth',2);
                    ylim([-1.5 1.5]);
%                     legend('Filtered Signal','First Derivative','Second Derivative'); legend box off;
                    xlabel('Time (sec)');
                    ylabel('Normalized Intensity');
                    set(gca,'FontSize',14);
                    set(gcf,'Position',[205,138,926,652]);
                    drawnow;
                    hold off; 
                end
            end
        end
%         totalNoise(pp) = sum(abs(rawSignal./max(rawSignal(:)) - tmpData));
    end
end

%%

riseIntDiff = riseEndInt-riseStrtInt;
fallIntDiff = fallStrtInt - fallEndInt;
tsRiseCleaned = tsRiseOT2_1;
tsFallCleaned = tsFallOT2_2;
sigWidthCleaned = sigWidth;
ampsCleaned = sigAmp;
onTransCleaned = sigsToMeasure;
onTransCleanedRaw = rawSigsToMeasure;
cellIDsCleaned = cellIDs;

% Remove points that don't have rising signal
tsRiseCleaned(riseIntDiff<=0) = [];
tsFallCleaned(riseIntDiff<=0) = [];
sigWidthCleaned(riseIntDiff<=0) = [];
ampsCleaned(riseIntDiff<=0) = [];
cellIDsCleaned(riseIntDiff<=0) = [];
onTransCleaned(riseIntDiff<=0,:) = [];
onTransCleanedRaw(riseIntDiff<=0,:) = [];

% Remove points that don't have falling signal
fallIntDiff(riseIntDiff<=0) = [];
tsRiseCleaned(fallIntDiff<=0) = [];
tsFallCleaned(fallIntDiff<=0) = [];
sigWidthCleaned(fallIntDiff<=0) = [];
ampsCleaned(fallIntDiff<=0) = [];
cellIDsCleaned(fallIntDiff<=0) = [];
onTransCleaned(fallIntDiff<=0,:) = [];
onTransCleanedRaw(fallIntDiff<=0,:) = [];

% Remove points that weren't analyzed
tsRiseCleaned(tsRiseCleaned<=0) = [];
sigWidthCleaned(tsFallCleaned<=0) = [];
ampsCleaned(tsFallCleaned<=0) = [];
cellIDsCleaned(tsFallCleaned<=0) = [];
onTransCleaned(tsFallCleaned<=0,:) = [];
onTransCleanedRaw(tsFallCleaned<=0,:) = [];
tsFallCleaned(tsFallCleaned<=0) = [];

cellIDsCleaned = ['Time' cellIDsCleaned]; 

sigStats = [mean(tsRiseCleaned) std(tsRiseCleaned); mean(tsFallCleaned) std(tsFallCleaned); mean(sigWidthCleaned) std(sigWidthCleaned);
    mean(ampsCleaned) std(ampsCleaned)];

% Display Histograms of Rise/Fall Half-Time
figure; histogram(tsRiseCleaned,'binwidth',0.5); hold on; histogram(tsFallCleaned,'binwidth',0.5); 
    xlabel('Time (sec)');
    ylabel('Frequency');
    title('\tau_{1/2} - ON-Transient');
    legend('Rising','Falling'); legend boxoff;
    set(gca,'FontSize',14,'fontweight','bold');
    if saveFlag
        saveas(gcf, [saveDir '/responseTimes.tiff']);
    end
    
% Display Signal Width Histogram
figure; histogram(sigWidthCleaned,'binwidth',0.5);
    xlabel('Time (sec)');
    ylabel('Frequency');
    title('Signal Width - ON-Transient');
    set(gca,'FontSize',14,'fontweight','bold');
    if saveFlag
        saveas(gcf, [saveDir '/peakWidth.tiff']);
    end
    
%Display Signal Amplitude Histogram
figure; histogram(ampsCleaned);
    xlabel('Amplitude');
    ylabel('Frequency');
    title('Signal Amplitude - ON-Transient');
    set(gca,'FontSize',14,'fontweight','bold');
    if saveFlag
        saveas(gcf, [saveDir '/peakAmplitude.tiff']);
    end
    
%% Save Quick-Look Stats

fid = fopen([saveDir '\results.txt'],'w');
fprintf(fid,'Cell Type: ON-Transient\n');
fprintf(fid,'N = %d\n', length(tsFallCleaned));
fprintf(fid,'t_rise = %0.2f +/- %0.2f seconds\n',sigStats(1,1),sigStats(1,2));
fprintf(fid,'t_fall = %0.2f +/- %0.2f seconds\n',sigStats(2,1),sigStats(2,2));
fprintf(fid,'sigWidth = %0.2f +/- %0.2f seconds\n',sigStats(3,1),sigStats(3,2));
fprintf(fid,'sigAmp = %0.2f +/- %0.2f\n',sigStats(4,1),sigStats(4,2));
fclose(fid);

%% Save Inidividual Cell Data

T = num2cell([tvals'; onTransCleanedRaw]');
T = table([cellIDsCleaned; T]);
tmpFname = 'rawAnalzyedSignals.csv';
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);

cellIDsCleaned{1} = 'Stats';
T = num2cell([tsRiseCleaned; tsFallCleaned; sigWidthCleaned; ampsCleaned]);
T = table([cellIDsCleaned; [{'Rise Time (sec)' 'Fall Time (sec)' 'Signal Width (sec)' 'Signal Amplitude'}' T]]);
tmpFname = 'analyzedSignalStats.csv';
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);