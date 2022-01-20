function [] = off_trans_Analysis(sigsToMeasure, rawSigsToMeasure, tvals, cellIDs, approxStart, approxEnd, viewPlots, saveFlag, saveDir)

%% Initialize varibles
[~, approxStartInd] = min(abs(tvals-approxStart)); % Index value corresponding to start of stimulus
[~, approxEndInd] = min(abs(tvals-approxEnd)); % Index value corresponding to end of stimulus


tsRiseOff = zeros(1,size(sigsToMeasure,1)); % rise times
tsFallOff = zeros(1,size(sigsToMeasure,1)); % fall times
sigStart1 = zeros(1,size(sigsToMeasure,1)); % start of rise signal
sigEnd1 = zeros(1,size(sigsToMeasure,1)); % end of rise signal
sigStart2 = zeros(1,size(sigsToMeasure,1)); % start of fall signal
sigEnd2 = zeros(1,size(sigsToMeasure,1)); % end of fall signal
sigWidth = zeros(1,size(sigsToMeasure,1)); % width of signal

%% Extract measurements from each signal

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
    
    if ~isempty(dfPk2)
        filtSigs(:,pp) = tmpData;
        % locate troughs in first derivative 
        troughData = islocalmin(dfData);
        troughInds = 1:length(troughData);
        troughInds = troughInds(troughData);
        zx = [zx; troughInds'];
        zx = unique(zx);

        sigStartG = (dfPk2 - dfTimes(zx))';
        sigStartG(sigStartG<0) = inf;
        [~,startInd] = min(sigStartG);

        % Calculate beginning and end of rise
        sigEnd1(pp) = dfTimes(zx(startInd+1));
        sigStart1(pp) = approxEnd;
        
        if approxEndInd < (zx(startInd+1)+1)
        
            [halfRiseTime, riseStrtInt(pp), riseEndInt(pp)] = calcHalfRise(tmpData(approxEndInd:zx(startInd+1)+1),tvals(approxEndInd:zx(startInd+1)+1));
            tsRiseOff(pp) = halfRiseTime - approxEnd;

            peakData = islocalmin(-dfData);
            peakInds = 1:length(peakData);
            peakInds = peakInds(peakData);
            zx2 = [zx; peakInds'];
            zx2 = unique(zx2);

            sigEndG = (sigEnd1(pp) - tvals(sigTroughInds));
            sigEndG(sigEndG>=0) = inf;
            [~,sigEndInd] = min(abs(sigEndG));

            % Calculate beginning and end of rise
            sigEnd2(pp) = tvals(sigTroughInds(sigEndInd));
            sigStart2(pp) = sigEnd1(pp);
%             tsFallOff_2(pp) = (sigEnd2(pp)-sigStart2(pp))/2;
            
            if (zx(startInd+1)+1) < sigTroughInds(sigEndInd)
                [halfFallTime, fallStrtInt(pp), fallEndInt(pp)] = calcHalfDecay(tmpData(zx(startInd+1)+1:sigTroughInds(sigEndInd)),tvals(zx(startInd+1)+1:sigTroughInds(sigEndInd)));
                tsFallOff(pp) = halfFallTime - sigEnd1(pp);
                
                sigWidth(pp) = halfFallTime - halfRiseTime;

                if viewPlots
                    plot(tvals,tmpData,'LineWidth',2); % plot signal
                    hold on;
%                     plot(tvals, rawSignal./max(abs(rawSignal(:))));
%                     plot(dfTimes,dfData./max(abs(dfData(:))),'LineWidth',2); % plot first derivative
%                     plot(df2Times,df2Data./max(abs(df2Data(:))),'LineWidth',2); 
                    plot([sigStart1(pp) sigStart1(pp)], [-1.5 1.5],'--k','LineWidth',2);
                    plot([sigStart2(pp) sigStart2(pp)], [-1.5 1.5],'--k','LineWidth',2);
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
    end
end

%% Display Histograms of Rise Half-Time

riseIntDiff = riseEndInt - riseStrtInt;
fallIntDiff = fallStrtInt - fallEndInt;
tsRiseCleaned = tsRiseOff;
tsFallCleaned = tsFallOff;
sigWidthCleaned = sigWidth;
offCleaned = sigsToMeasure;
offCleanedRaw = rawSigsToMeasure;
cellIDsCleaned = cellIDs;

% Remove points that don't have rising signal
tsRiseCleaned(riseIntDiff<=0) = [];
tsFallCleaned(riseIntDiff<=0) = [];
sigWidthCleaned(riseIntDiff<=0) = [];
cellIDsCleaned(riseIntDiff<=0) = [];
offCleaned(riseIntDiff<=0,:) = [];
offCleanedRaw(riseIntDiff<=0,:) = [];

% Remove points that don't have falling signal
fallIntDiff(riseIntDiff<=0) = [];
tsRiseCleaned(fallIntDiff<=0) = [];
tsFallCleaned(fallIntDiff<=0) = [];
sigWidthCleaned(fallIntDiff<=0) = [];
cellIDsCleaned(fallIntDiff<=0) = [];
offCleaned(fallIntDiff<=0,:) = [];
offCleanedRaw(fallIntDiff<=0,:) = [];

% Remove points that weren't analyzed
tsRiseCleaned(tsRiseCleaned<=0) = [];
sigWidthCleaned(tsFallCleaned<=0) = [];
cellIDsCleaned(tsFallCleaned<=0) = [];
offCleaned(tsFallCleaned<=0,:) = [];
offCleanedRaw(tsFallCleaned<=0,:) = [];
tsFallCleaned(tsFallCleaned<=0) = [];

cellIDsCleaned = ['Time' cellIDsCleaned]; 

sigStats = [mean(tsRiseCleaned) std(tsRiseCleaned); mean(tsFallCleaned) std(tsFallCleaned); mean(sigWidthCleaned) std(sigWidthCleaned)];

% Display Histograms of Rise/Fall Half-Time
figure; histogram(tsRiseCleaned,'binwidth',1); hold on; histogram(tsFallCleaned,'binwidth',1); 
    xlabel('Time (sec)');
    ylabel('Frequency');
    title('OFF \tau_{1/2}');
    legend('Rising','Falling'); legend boxoff;
    set(gca,'fontweight','bold','FontSize',14);
    if saveFlag
        saveas(gcf, [saveDir '/responseTimes.tiff']);
    end
    
% Display Signal Width Histogram
figure; histogram(sigWidthCleaned,'binwidth',1);
    xlabel('Time (sec)');
    ylabel('Frequency');
    title('OFF Signal Width');
    set(gca,'fontweight','bold','FontSize',14);
    if saveFlag
        saveas(gcf, [saveDir '/peakWidth.tiff']);
    end
    
%% Save quick look text file
fid = fopen([saveDir '\results.txt'],'w');
fprintf(fid,'Cell Type: OFF\n');
fprintf(fid,'N = %d\n', length(offCleaned));
fprintf(fid,'t_rise = %0.2f +/- %0.2f seconds\n',sigStats(1,1),sigStats(1,2));
fprintf(fid,'t_fall = %0.2f +/- %0.2f seconds\n',sigStats(2,1),sigStats(2,2));
fprintf(fid,'sigWidth = %0.2f +/- %0.2f seconds\n',sigStats(3,1),sigStats(3,2));
fclose(fid);

%% Save Inidividual Cell Data
T = num2cell([tvals'; offCleanedRaw]');
T = table([cellIDsCleaned; T]);
tmpFname = 'rawAnalzyedSignals.csv';
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);

cellIDsCleaned{1} = 'Stats';
T = num2cell([tsRiseCleaned; tsFallCleaned; sigWidthCleaned]);
T = table([cellIDsCleaned; [{'tRise (sec)' 'tFall (sec)' 'sigWidth (sec)'}' T]]);
tmpFname = 'analyzedSignalStats.csv';
writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);