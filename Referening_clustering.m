warning ('off','all')
ClusterName=["ON-Sust1","ON-Sust2","ON-Sust3","ON-Trans", "ON-OFF","OFF-Trans","OFF-Sust","OFF-Supp1","OFF-Supp2","NR","OFF-TRANS2","ON-OFF-TRANS"];

%F = uigetdir('','Select Input-folder');
F='H:\Liang Li data\Ca imaging of RGC and analysis by GCaMP\Analysis test and optimization\!!!Shaobo Data Analysis\NNRI\10. regroup standard of Naive retina  !!!';
%Above Dir path is the grouping result of the naive RGC, which is for
%reference of the denegerative RGC grouping.
FileList = dir(fullfile(F, '**', '*.csv'));
for iFile = 1:numel(FileList)
    file = fullfile(FileList(iFile).folder, FileList(iFile).name);
    FItenAll = csvread(file);
    FItenAll=FItenAll(:,2:end);
    MeanIn = mean(FItenAll');
    
    MeanIn = MeanIn/norm(MeanIn);
    MeanInten(:,iFile)=MeanIn;

    %figure,
    subplot(2,5,iFile);
    hold on;
    plot(FItenAll,'g');
    plot(MeanInten(:,iFile));
    title(['Cluster' num2str(iFile) ClusterName(iFile)]);
    %figure, plot(MeanInten(:,iFile));
end

[file,path]  = uigetfile('*.csv','where is the csv file of the traces?');
FItenAll = csvread(fullfile(path,file));
FItenAll = FItenAll-1;

cutoffs = input('Do you want to filter out the traces with less than 0.15~-0.20 variation?  (yes) or (no):  ','s');
pcutoff = 0.15;
ncutoff = -0.20;  
% positive cutoff is 0.15, changable; negative cutoff 0.2 changable.
if cutoffs(1) == 'y'
    disp('Removing the cells with less than 20% variation')
    varia=[min(FItenAll,[],2),max(FItenAll,[],2)];
    %varia = max(varia,[],2);
    NoresList = find(varia(:,1) > ncutoff & varia(:,2) < pcutoff);
    Noresponse = FItenAll(NoresList,:);
    disp(['Removed  ' num2str(size(NoresList,1)) '  Traces']);
    % Noresponse correspond to the filtered traces (considered as No responding cells)
    FItenAll(NoresList,:)=[];
    clear varia;
    
    csvwrite('NoReponse-Trace-File.csv',Noresponse');
    csvwrite('Response-Trace-File.csv',FItenAll');
    
end

% FItenAllRaw=FItenAllRaw';
%%

FItenAll=FItenAll';
figure
Threshold = 0.6; % threshold for correlation
CellID=[];

for i=1:size(FItenAll,2)

    Trace=FItenAll(:,i);
    Trace = Trace/norm(Trace);
    tmpData = lowpass(Trace,1e-10);

    FilteredData(:,i)=tmpData;


    for j=1:size(MeanInten,2)
        R = xcorr(FilteredData(:,i),MeanInten(:,j),0);
        RR(j)= R;
    end
    [M,ClusterID] = max(RR);
    plot(Trace);

    %     tmpData = tmpData./max(abs(tmpData(:)));

    hold on
    plot(tmpData,'r');
    if M>Threshold
        title(['Cell' num2str(i) ' ' ClusterName(ClusterID)  ' Correlation Ratio: ' num2str(M)]);
        plot(MeanInten(:,ClusterID),'g');
    else
        title(['Cell' num2str(i) ' UnAssigned(Probably' ClusterName(ClusterID)  ') Correlation Ratio: ' num2str(M)]);
        ClusterID = 10;
    end
    hold off
    pause(0.001);
    CellID(i)=ClusterID;
end

% CellID(find(CellID==11))=6;
% CellID(find(CellID==12))=5;
%(fO)

% for i=1:size(FItenAll,2)
%     
%     Trace=FItenAll(:,i);
%     for j=1:size(MeanInten,2)
%         R = corrcoef(Trace,MeanInten(:,j));
%         RR(j)= R(1,2);
%     end
%     [M,ClusterID] = max(RR);
%     plot(Trace);
%     
%     
%     %     tmpData = tmpData./max(abs(tmpData(:)));
%     tmpData = lowpass(Trace,1e-10);
%     FilteredData(:,i)=tmpData;
%     hold on
%     plot(tmpData,'r');
%     if M>Threshold
%         title([ClusterName(ClusterID)  ' Correlation Ratio: ' num2str(M)]);
%     else
%         
%         title(['UnAssigned(Probably' ClusterName(ClusterID)  ') Correlation Ratio: ' num2str(M)]);
%         ClusterID = 11;
%     end
%     hold off
%     pause(0.1); % display time adjustable 'pause(0.5)'
%     CellID(i)=ClusterID;
% end


%%
viewPlots = 1;
saveFlag = 1;
approxStart = 10;
approxEnd = 30;
tvals = linspace(0,60,283);
tvals = tvals(4:279);
tvals=tvals';

sigsToMeasure = zeros(1,size(FItenAll,2));
rawSigsToMeasure = zeros(1,size(FItenAll,2));
cellIDsToMeasure = {};

sigsToMeasure(1,:) = [];
rawSigsToMeasure(1,:) = [];


CellI=CellID;

for kk = [1:10] %Can be changed for specific cluster by typing in related cluster number
    %     try
    
    if size(find(CellI==kk),2)>0
        figure
        plot(FItenAll(:,CellI==kk),'g');
        xticks([1,45,92,139,186,233,280]);
        xticklabels({'0','10','20','30','40','50','60'});
        %axis([0,283,Ymin,Ymax]);
        MeanInten(:,kk)=mean(FItenAll(:,CellI==kk)');
        hold on
        plot(MeanInten(:,kk),'r');
        title(['Group' num2str(kk)])
        hold off;
        
        cellType = input('Please Specify Cell Type: ','s');
        %kk;
        tmpSigs = FItenAll(:,CellI==kk);
        %sigsToMeasure = tmpSigs;
        %         disp('Filtering Data');
        %         for pp = 1:size(tmpSigs,1)
        %             % denoise data
        %             tmpData = tmpSigs(pp,:);
        %             %     tmpData = tmpData./max(abs(tmpData(:)));
        %             tmpData = lowpass(tmpData,1e-10);
        %             sigsToMeasure(pp,:) = tmpData;
        %         end
        %         disp('Filtering Done!');
        tmpSigsRaw = FItenAll(:,CellI==kk);
        
        rawSigsToMeasure = tmpSigsRaw;
        %figure,plot(rawSigsToMeasure);
        cellIDsToMeasure = num2cell(find(CellI==kk));
        %cellIDsToMeasure=cellIDsToMeasure';
        sigsToMeasure = FItenAll(:,CellI==kk);
        %figure,plot(sigsToMeasure);
        % %         sigsToMeasure = FilteredData;
       disp(['Cell Number: ' num2str(size(find(CellI==kk),2)) ]);
        
        if saveFlag
            svCheck = 1;
            while exist([path file(1:end-4) '_' 'Analysis_' num2str(svCheck)],'dir')
                svCheck = 1 + svCheck;
            end
            foldername = ['Cluster' num2str(kk)];
            
            mkdir([path file(1:end-4) '_' 'Analysis_' foldername]);
            
            saveDir = [path file(1:end-4) '_' 'Analysis_' foldername];
        end
        
        
        
        if strcmp(cellType,'OFF-TRANS')
            off_trans_Analysis(sigsToMeasure',rawSigsToMeasure',tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
        elseif strcmp(cellType, 'ON-SUST1')
            on_sust1_Analysis(sigsToMeasure',rawSigsToMeasure',tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
        elseif strcmp(cellType, 'ON-SUST2')
            on_sust2_Analysis(sigsToMeasure',rawSigsToMeasure',tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
        elseif strcmp(cellType, 'ON-SUST3')
            on_sust3_Analysis(sigsToMeasure',rawSigsToMeasure',tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
        elseif strcmp(cellType, 'ON-OFF')
            on_off_Analysis(sigsToMeasure',rawSigsToMeasure',tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
        elseif strcmp(cellType, 'ON-TRANS')
            figure;
            on_trans_Analysis(sigsToMeasure',rawSigsToMeasure',tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
        elseif strcmp(cellType, 'SUPP2')
            supp2_Analysis(sigsToMeasure',rawSigsToMeasure',tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
        elseif strcmp(cellType, 'OFF-SUST')
            off_sust_Analysis(sigsToMeasure',rawSigsToMeasure',tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
        elseif strcmp(cellType, 'SUPP1')
            supp1_Analysis(sigsToMeasure',rawSigsToMeasure',tvals,cellIDsToMeasure,approxStart,approxEnd,viewPlots,saveFlag,saveDir);
        else
            disp('Cell Type Not Recognized');
            tmpFname='Unassigned.csv';
            csvwrite([saveDir '\' tmpFname],rawSigsToMeasure);
        end
    end
    %     catch
    %         disp(['Error found in Cluster ' num2str(kk)]);
    %     end
end