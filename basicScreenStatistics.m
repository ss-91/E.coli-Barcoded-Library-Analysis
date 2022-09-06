%% This script plots basic statistics of screen results (no biological anaysis)
%% user definitions
fileToSave = 'RPM_COUNTS_QS10.mat'; % name of output file
fileToLoad = 'dataset_QS10.mat'; % this file is the output of fastq2barcodeCounts.m
load dataset_QS10; 

% Two options for specifying data structure in fileToLoad  (conditions, replicates) 
% Option 1 (index automatically by number of conditions and replicates)
nReplicates = 3;
nConditions = 2;
replicateInx = zeros(nConditions,nReplicates)'; % replicate index appear in each column. This can also be added manually
replicateInx(:) = [1:nReplicates*nConditions]; replicateInx = replicateInx'; % order 1:n number by replicate/condition order

% Option 2 (manually specify the indexing)
%replicateInx = [5,6;1,2;3,4;11,12;7,8;9,10]; % OR - specify the indexes (replicates in the same row, conditions in separate rows)
%[nConditions,nReplicates] = size(replicateInx);

nScreens = length(replicateInx(:));


%% make sensible names
for i=1:nScreens % should match size of dataset
    curFileName = dataset(i).fileName;
    newStr = split(curFileName,["_",'.']); % split by _ or .
    datasetName{i} = newStr{1};
end

%% calucatle RPM, COUNTS and compactRPM
RPM = []; % master matrix to hold all RPMs
compactRPM =[];
COUNTS = [];
labels = {}; % labels for RPM matrix
compactLabels = {};

for iCond=1:nConditions
    for iRep=1:nReplicates
        iDataset = replicateInx(iCond,iRep);
        RPM = [RPM,dataset(iDataset).hits.RPM];
        COUNTS = [COUNTS,dataset(iDataset).hits.counts];
        labels = {labels{:},datasetName{iDataset}};
    end
end

% calculate compactRPM (RPM that averages over replicates)
for iCond=1:nConditions
    curSetInx = replicateInx(iCond,:);
    miniRPM = [];
    iDataset = replicateInx(iCond,iRep);
    for iRep=1:nReplicates
        miniRPM = [miniRPM,dataset(replicateInx(iCond,iRep)).hits.RPM];
    end
    compactRPM = [compactRPM,nanmean(miniRPM,2)];
    compactLabels = {compactLabels{:},datasetName{iDataset}}; % label by the last replicate of this condition
end

%% save everything

save(fileToSave,'RPM','compactRPM','COUNTS','labels','compactLabels');

%% Now for some plot

%% plot the number of identified/unidentified barcode per index
figure; hold on;
subplot(4,1,[1:3]);
%inx = replicateInx(:);
inx=[1,2,3,4,5,6]
bar([[dataset(inx).nPositiveBarcodes];[dataset(inx).nNegativeBarcodes]]','stacked');
set(gca,'xtick',[1:length(inx)],'xticklabel',datasetName(inx));
xtickangle(45); grid on;
ylabel('Number of reads');
legend({'Reads with barcode','Reads w/o barcodes'});
title('Number of identified barcodes');
subplot(4,1,4);
r = [dataset(inx).nPositiveBarcodes]./[[dataset(inx).nPositiveBarcodes] + [dataset(inx).nNegativeBarcodes]];
bar(r,'facecolor','k');
set(gca,'xtick',[1:length(inx)],'xticklabel',datasetName(inx),'ylim',[0 1]);
xtickangle(45); grid on;
ylabel('pos reads / total reads');

%% plot replicates as scatters

figure; hold on;
i = 0;
for iCond=1:nConditions
    for iRep = 1:nReplicates
        for jRep = (iRep+1):nReplicates
            i = i+1;
            subplot(nConditions,nchoosek(nReplicates,2),i); hold on;
            inxDataset = replicateInx(iCond,:);
            box on;
            plot(log10(dataset(inxDataset(iRep)).hits.RPM),log10(dataset(inxDataset(jRep)).hits.RPM),'.k');
            set(gca,'xlim',[0 5],'ylim',[0 5]);
            axis square; grid on;
            xlabel(datasetName(inxDataset(iRep))); ylabel(datasetName(inxDataset(jRep)));
            % building the RPM matrix as we go
        end
    end
end

suptitle('Biological repeats - RPM');


%% plot correlations between individual screens (replicates and conditions)
y = log2(RPM);
y(RPM==0)=nan;

figure;
subplot(1,2,1);
[r p] = corr(y,'Type','Pearson','Rows','complete');
imagesc(p); colormap(jet);
set(gca,'xtick',[1:length(labels)],'ytick',[1:length(labels)],'xticklabel',labels,'yticklabel',labels);
xtickangle(90);
axis square;
title('Pearson corr (linear)');

subplot(1,2,2);
[r p] = corr(y,'Type','Spearman','Rows','complete');
imagesc(p); colormap(gray);
set(gca,'xtick',[1:length(labels)],'ytick',[1:length(labels)],'xticklabel',labels,'yticklabel',labels);
xtickangle(90);
axis square;
title('Spearman corr (rank)');

%% Cluster by correlations and plot as clustergram
[r p] = corr(y,'Type','Pearson','Rows','complete');
cg = clustergram(r,'rowlabels',labels,'columnlabels',labels,'linkage','complete','Symmetric',false,'RowPDist','euclidean','ColumnPDist','euclidean')
cg.Colormap = gray;
cg.DisplayRange = 2;
cg.DisplayRatio = [0.1 0.1];
cgAxes =plot(cg); % Use the plot function to plot to a separate figure and output the axes
set(cgAxes, 'Clim', [0.5,1]) % Set colour limit or other axes properties.





