%% Initialize
clear; clc; close all;

% dir0 ='G:\Liang Li data\Ca imaging of RGC and analysis by GCaMP\Analysis test and optimization\!Colaborator David Miller analysis\RGC Analysis Software v3\';
% addpath(dir0)

viewPlots = 1;
saveFlag = 1;
approxStart = 10;
approxEnd = 30;
%% Generate Colormap Values

cmapL = 256;

c1 = [255 0 0]/255;
blck = [0 0 0];
c12Black = zeros(cmapL/2,1,3);
c12Black(:,1,:) = [linspace(c1(1),blck(1),cmapL/2)' linspace(c1(2),blck(2),cmapL/2)' linspace(c1(3),blck(3),cmapL/2)'];

c2 = [20 255 0]/255;
c22Black = zeros(cmapL/2,1,3);
c22Black(:,1,:) = [linspace(c2(1),blck(1),size(c22Black,1))' linspace(c2(2),blck(2),size(c22Black,1))' linspace(c2(3),blck(3),size(c22Black,1))'];

cmap = [c12Black; flipud(c22Black)];
cmap = squeeze(cmap);
%% Open File

[fname, dir0] = uigetfile('*.csv');
dataIn = readtable([dir0 fname],'ReadVariableNames',0);
cellIDs = table2array(dataIn(1,3:end));
tvals = table2array(dataIn(2:end,2));
dataIn = table2array(dataIn(2:end,3:end));
tvals = str2double(tvals);
dataIn = str2double(dataIn);

dataIn = dataIn';
filtSig = dataIn;
[~, approxStartInd] = min(abs(tvals-approxStart));
[~, approxEndInd] = min(abs(tvals-approxEnd));
%% Make Save Directory

clc; 
if saveFlag
    svCheck = 1;
    while exist([dir0 fname(1:end-4) '_' 'Analysis_' num2str(svCheck)],'dir')
        svCheck = 1 + svCheck;
    end
    
    mkdir([dir0 fname(1:end-4) '_' 'Analysis_' num2str(svCheck)]);
    
    saveDir = [dir0 fname(1:end-4) '_' 'Analysis_' num2str(svCheck)];
end
%% Denoise data

clc;
disp('Denoising Data (may take some time)');
for pp = 1:size(dataIn,1)
     % denoise data
    tmpData = dataIn(pp,:);
%     tmpData = tmpData./max(abs(tmpData(:)));
    tmpData = lowpass(tmpData,1e-10); 
    filtSig(pp,:) = tmpData;
end
%% Display Final Signals

figure(1); imagesc(tvals,1:size(filtSig,1),filtSig./max(abs(filtSig),[],2)); 
        colormap(cmap); caxis([-1 1]); colorbar; 
        xlabel('Time (sec)');
        ylabel('Signal #');
        set(gca,'fontweight','bold','FontSize',14);
        set(gcf,'Position',[680,61,560,917]);
        hold on; drawnow;
        
saveas(gcf,[saveDir '\' 'original_Signals.tiff']);
%% Cluter Signals

clc;
clustNo = input('How many clusters?: ');
s = mdwtcluster(filtSig(:,approxStartInd:end-20)./max(abs(filtSig),[],2),'maxClust',clustNo);
%% Reject Misclassified Signals

idxCLU = s.IdxCLU;
reclusteredSignals = zeros(1,size(filtSig,2));
reclusteredCellIDs = {};
clustBounds = [];

for pp = 1:clustNo
    tmpSigs = filtSig(idxCLU(:,1)==pp,:);
    reclusteredSignals(end+1:end+size(tmpSigs,1),:) = tmpSigs;
    clustBounds(end+1) = size(reclusteredSignals,1)-1;
    tmpIDs = cellIDs(idxCLU(:,1)==pp);
    reclusteredCellIDs(end+1:end+size(tmpSigs,1)) = tmpIDs;
end
reclusteredSignals(1,:) = [];

figure(6); 

subplot(121); imagesc(tvals,1:size(reclusteredSignals,1),reclusteredSignals./max(abs(reclusteredSignals),[],2)); 
        colormap(cmap); caxis([-1 1]); colorbar; 
        xlabel('Time (sec)');
        ylabel('Signal #');
        set(gca,'fontweight','bold','FontSize',14);
        title('Clustered Signals'); hold on; drawnow;

for pp = 1:length(clustBounds)-1
    plot([tvals(1) tvals(end)], [clustBounds(pp) clustBounds(pp)], ':w','Linewidth',2);
end

for pp = 1:clustNo
    if pp == 1
        text(tvals(25),(clustBounds(pp)+1)/2,num2str(pp),'Color','White');
    elseif pp==clustNo
        text(tvals(25),(size(reclusteredSignals,1)+clustBounds(pp-1))/2,num2str(pp),'Color','White');
    else
        text(tvals(25),(clustBounds(pp)+clustBounds(pp-1))/2,num2str(pp),'Color','White');
    end
end
plot([approxStart approxStart], [1 size(reclusteredSignals,1)],'--y');
plot([approxEnd approxEnd], [1 size(reclusteredSignals,1)],'--y');
set(gca,'fontweight','bold','FontSize',14,'YDir','Normal');
hold off;

subplot(122); hold on; 
clustSignals = zeros(clustNo,size(reclusteredSignals,2));
for pp = 1:clustNo
    tmpSigs = filtSig(idxCLU(:,1)==pp,:);
    clustSignals(pp,:) = mean(tmpSigs)./(max(abs(mean(tmpSigs)))*2);
    plot(tvals, clustSignals(pp,:)+pp,'LineWidth',2);
end

for aa = 1:size(clustSignals,1)
    sigCent = aa;
    plot([tvals(1) tvals(end)], [sigCent sigCent], '--k');
end
plot([approxStart approxStart], [0 (size(clustSignals,1)+0.5)],'--r');
plot([approxEnd approxEnd], [0 (size(clustSignals,1)+0.5)],'--r');
ylim([0 size(clustSignals,1)+1]);
xlabel('Time (sec)');
ylabel('Cluster #');
title('Mean Cluster Traces');
set(gca,'YTick',1:size(clustSignals,1),'fontweight','bold','FontSize',14);
set(gcf,'Position',[680,61,1031,917]);
saveas(gcf,[saveDir '\' 'clustered_Signals.tiff']);
%% Save Traces

for pp = 1:clustNo
    rawClustData = dataIn(idxCLU(:,1)==pp,:);
    tmpIDs = reclusteredCellIDs(idxCLU(:,1)==pp);
    tmpIDs = ['_' 'Time' tmpIDs];
    T = table([tmpIDs; num2cell([(1:length(tvals)); tvals'; rawClustData])']);
    tmpFname = ['cluster_' num2str(pp,'%02.f') '.csv'];
    writetable(T,[saveDir '\' tmpFname],'WriteVariableNames',0);
end
%% Ask user to input which signals to analyze

clc; 
kSigs = input('Enter signals to keep (use brackets): ');
%% Extract Signals to analyze

sigsToMeasure = zeros(1,size(reclusteredSignals,2));
rawSigsToMeasure = zeros(1,size(reclusteredSignals,2));
cellIDsToMeasure = {};
for kk = kSigs
    tmpSigs = filtSig(idxCLU(:,1)==kk,:);
    sigsToMeasure(end+1:end+size(tmpSigs,1),:) = tmpSigs;
    tmpSigsRaw = dataIn(idxCLU(:,1)==kk,:);
    rawSigsToMeasure(end+1:end+size(tmpSigs,1),:) = tmpSigsRaw;
    tmpIDs = reclusteredCellIDs(idxCLU(:,1)==kk);
    cellIDsToMeasure(end+1:end+length(tmpIDs)) = tmpIDs;
end
sigsToMeasure(1,:) = [];
rawSigsToMeasure(1,:) = [];

figure(7); imagesc(tvals,1:size(sigsToMeasure,1),(sigsToMeasure./max(abs(sigsToMeasure),[],2))); colormap(cmap); caxis([-1 1]); colorbar; 
    xlabel('Time (sec)'); ylabel('Signal #'); title('Reclassified Signals');
    set(gca,'fontweight','bold', 'FontSize',14); set(gcf,'Position',[680,61,560,917]); drawnow;
    
saveas(gcf,[saveDir '\' 'Analyzed_Signal_Traces.tiff']);
%% Display Cluster Traces

figure; hold on;
gplot = 0;
rplot = 0;
clustSignals = zeros(clustNo,size(reclusteredSignals,2));
for pp = 1:clustNo
    tmpSigs = filtSig(idxCLU(:,1)==pp,:);
    if sum(pp==kSigs)
        gplot = gplot + 1;
        clustSignals(pp,:) = mean(tmpSigs)./(max(abs(mean(tmpSigs)))*2);
        if gplot == 1
            ln1 = plot(tvals, clustSignals(pp,:)+pp,'g','LineWidth',2);
        else
            plot(tvals, clustSignals(pp,:)+pp,'g','LineWidth',2);
        end
    else
        rplot = rplot + 1;
        clustSignals(pp,:) = mean(tmpSigs)./(max(abs(mean(tmpSigs)))*2);
        if rplot == 1
            ln2 = plot(tvals, clustSignals(pp,:)+pp,'r','LineWidth',2);
        else
            plot(tvals, clustSignals(pp,:)+pp,'r','LineWidth',2);
        end
    end
end

for aa = 1:size(clustSignals,1)
    sigCent = aa;
    plot([tvals(1) tvals(end)], [sigCent sigCent], '--k');
end
plot([approxStart approxStart], [0 (size(clustSignals,1)+0.5)],'--r');
plot([approxEnd approxEnd], [0 (size(clustSignals,1)+0.5)],'--r');
ylim([0 size(clustSignals,1)+1]);
xlabel('Time (sec)');
ylabel('Cluster #');
set(gca,'YTick',1:size(clustSignals,1),'fontweight','bold','FontSize',14);
set(gcf,'Position',[680,61,560,917]);    
if rplot == 0 && gplot > 0
    legend(ln1,{'Analyzed Signals'},'Location','Northwest'); legend boxoff;
elseif rplot > 0 && gplot == 0
    legend(ln2,{'Removed Signals'},'Location','Northwest'); legend boxoff;
else
    legend([ln1 ln2],{'Analyzed Signals','Removed Signals'},'Location','Northwest'); legend boxoff;
end

saveas(gcf,[saveDir '\' 'Reclassified_Signal_Traces.tiff']);
%% Analyze Selected Signals
clc;
cellType = input('Please Specify Cell Type: ','s');

if strcmp(cellType,'OFF-TRANS')
    off_trans_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
elseif strcmp(cellType, 'ON-SUST1')
    on_sust1_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
elseif strcmp(cellType, 'ON-SUST2')
    on_sust2_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
elseif strcmp(cellType, 'ON-SUST3')
    on_sust3_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
elseif strcmp(cellType, 'ON-OFF')
    on_off_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
elseif strcmp(cellType, 'ON-TRANS')
    on_trans_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
elseif strcmp(cellType, 'SUPP2')
    supp2_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
elseif strcmp(cellType, 'OFF-SUST')
    off_sust_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
elseif strcmp(cellType, 'SUPP1')
    supp1_Analysis(sigsToMeasure,rawSigsToMeasure,tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
else
    disp('Cell Type Not Recognized');
end