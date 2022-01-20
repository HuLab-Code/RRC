clear; clc; close all;

dir0 = ['M:\Stanford\2020-12-15 Ca2+ activity data issues\Demo data and R code for heatmap\Heatmap of ONC and SOHU\Heatmaps\'];
addpath(dir0);

fileList = dir([dir0 '*.xlsx']);

for filNum = 1:size(fileList,1)
    close all;
    fname = fileList(filNum).name;
    cellData = xlsread([dir0 fname]);
    
    tvals = cellData(:,1);

    cellSigs = rot90(cellData(:,3:end));
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

    figure; imagesc(tvals,1:size(cellSigs,1),cellSigs./max(abs(cellSigs),[],2)); colormap(cmap); colorbar; caxis([-1 1]);
    set(gcf, 'Position', [680,54,560,924]);
    xlabel('Time (sec)');
    ylabel('Signal Number');
    set(gca,'FontSize',14);
%     drawnow;
%     saveas(gcf,[dir0 fname(1:end-4) '_unfilt.tiff']);

    %%
    approxStart = 10;
    approxEnd = 30;

    [~, approxStartInd] = min(abs(tvals-approxStart));
    [~, approxEndInd] = min(abs(tvals-approxEnd));
    onSust2 = cellSigs;

    filtSig = onSust2;

    for pp = 1:size(onSust2,1)
         % denoise data
        tmpData = onSust2(pp,:);
        tmpData = tmpData./max(abs(tmpData(:)));
        tmpData = lowpass(tmpData,1e-10); 
        filtSig(pp,:) = tmpData;
    end

    % figure; imagesc(tvals,1:size(onSust2,1),onSust2./max(abs(onSust2),[],2)); colormap(cmap); caxis([-1 1]);
    figure; imagesc(tvals,1:size(filtSig,1),filtSig./max(abs(filtSig),[],2)); 
        colormap(cmap); caxis([-1 1]); colorbar; 
        xlabel('Time (sec)');
        ylabel('Signal #');
        title('Original Signals');
        set(gca,'FontSize',14); 
        set(gcf, 'Position', [680,54,560,924]);
        drawnow;
    saveas(gcf,[dir0 fname(1:end-5) '_filt.tiff']);


    %% Save Imgs
    rawSigs = cellSigs./max(abs(cellSigs),[],2);
    filtSigs = filtSig./max(abs(filtSig),[],2);

    Numcolor = size(cmap, 1);
    Imgmin = min(rawSigs(:)) ;
    Imgmax = max(rawSigs(:)) ;
    mappedRawSigs = uint8((rawSigs-Imgmin)./(Imgmax-Imgmin).* (Numcolor-1) ) ;

    Imgmin = min(filtSigs(:)) ;
    Imgmax = max(filtSigs(:)) ;
    mappedFiltSigs = uint8((filtSigs-Imgmin)./(Imgmax-Imgmin).* (Numcolor-1) ) ;


%     imwrite(mappedRawSigs, (cmap), [dir0 fname(1:end-4) '_unfilt_noAx.tiff']);
    imwrite(mappedFiltSigs, (cmap), [dir0 fname(1:end-5) '_filt_noAx.tiff']);
end