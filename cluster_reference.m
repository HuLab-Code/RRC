warning ('off','all')
ClusterName=["ON-Sust2","OFF-Supp1","ON-Trans","ON-Sust1", "ON-OFF","OFF-Trans","OFF-Supp2","OFF-Sust","ON-Sust3","No Response"];

%F = uigetdir('','Select Input-folder');
F='H:\Liang Li data\Ca imaging of RGC and analysis by GCaMP\!!! Figures for paper\Figure1 Method\!!S-D method\2. cluster\10. regroup standard of Naive retina  !!!';
%Above Dir path is the grouping result of the naive RGC, which is for
%reference of the denegerative RGC grouping.
FileList = dir(fullfile(F, '**', '*.csv'));
figure
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
end

%%
FItenAllTrace=[];
F=uigetdir('','Select input folder');
FileList=dir(fullfile(F,'**','*.csv'));
for iFile = 1:numel(FileList)
    file=fullfile(FileList(iFile).folder,FileList(iFile).name);
    %[file,path]  = uigetfile('*.csv','where is the csv file of the traces?');
    %Table = readtable(fullfile(path,file));
    Table = readtable(file);
    cv = Table.Properties.VariableNames;
    FItenAllTraceT = table2array(Table);
    FItenAllTraceT = FItenAllTraceT(4:279,3:end);
    size(FItenAllTraceT)
    FItenAllTrace = [FItenAllTrace,FItenAllTraceT];
end
% FItenAll = FItenAll';
%FItenAll = FItenAll-1;

% cutoffs = input('Do you want to filter out the traces with less than 0.15~-0.20 variation?  (yes) or (no):  ','s');
% pcutoff = 0.15;
% ncutoff = -0.2;
% % positive cutoff is 0.1, changable; negative cutoff 0.2 changable.
% if cutoffs(1) == 'y'
%     disp('Removing the cells with less than 20% variation')
%     varia=[min(FItenAll,[],2),max(FItenAll,[],2)];
%     %varia = max(varia,[],2);
%     NoresList = find(varia(:,1) > ncutoff & varia(:,2) < pcutoff);
%     Noresponse = FItenAll(NoresList,:);
%     disp(['Removed  ' num2str(size(NoresList,1)) '  Traces']);
%     % Noresponse correspond to the filtered traces (considered as No responding cells)
%     FItenAll(NoresList,:)=[];
%     clear varia;
%
%     csvwrite('NoReponse-Trace-File.csv',Noresponse');
%     csvwrite('Response-Trace-File.csv',FItenAll');
%
% end

% FItenAllRaw=FItenAllRaw';
%%
%%
CellID=[];
Interval = input('How many time points have you measured for this experiment? (input like 5 or 6)   ');

FItenAll = FItenAllTrace(:,1:Interval:end);
%FItenAll=FItenAll';
figure
Threshold = 0.7; % threshold for correlation


for i=1:size(FItenAll,2)

    Trace=FItenAll(:,i);
    Trace = Trace/norm(Trace);
    tmpData = lowpass(Trace,1e-10);

    FilteredData(:,i)=tmpData;


    for j=1:size(MeanInten,2)
        R = xcorr(FilteredData(:,i),MeanInten(:,j),0,'normalized');
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
    pause; % display time adjustable 'pause(0.5)'
    CellID(i,1)=ClusterID;
end


FItenAll = FItenAllTrace(:,Interval:Interval:end);
%FItenAll=FItenAll';
figure
Threshold = 0.7; % threshold for correlation

for i=1:size(FItenAll,2)

    Trace=FItenAll(:,i);
    Trace = Trace/norm(Trace);
    tmpData = lowpass(Trace,1e-10);

    FilteredData(:,i)=tmpData;


    for j=1:size(MeanInten,2)
        R = xcorr(FilteredData(:,i),MeanInten(:,j),0,'normalized');
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
        ClusterID = 11;
    end
    hold off
    pause; % display time adjustable 'pause(0.5)'
    CellID(i,2)=ClusterID;
end

for i=1:10
    for j=1:10
        CellT= CellID(CellID(:,1)==i,:);
        data(i,j)=size(find(CellT(:,2)==j),1);
    end
end

figure
alluvialflow(data,ClusterName,ClusterName,'Cluster Change');