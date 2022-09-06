%% Plot pathway enrichment results

inputFileName = 'PathwayEnrichment.xlsx';
inxScreen = [1]; % a subset of the screen to use

%% load enrichment results (excel file)
[status,sheets] = xlsfinfo(inputFileName); % get the sheet names

for i=1:max(inxScreen) % iterate over sheets and get data
    [num,txt,raw] = xlsread(inputFileName,i);
    path{i} = txt(2:end,1);
    pval{i} = num(:,4);
    setSize{i} = num(:,5);
end

path = path(inxScreen); 
pval = pval(inxScreen);
setSize = setSize(inxScreen);
screenLabels = sheets(inxScreen);

uniquePath = unique(cat(1,path{:}));



%% split paths to kegg, go and ecocyc

pathGO = {};
pathKEGG = {};
pathEcoCyc = {};

for i=1:length(uniquePath)
    curPath = uniquePath{i};
    if(~isempty(regexp(curPath,'^GO:')))
        pathGO{end+1} = curPath;
    elseif(~isempty(regexp(curPath,'^eco')))
        pathKEGG{end+1} = curPath;
    else
        pathEcoCyc{end+1} = curPath;
    end
end

%% calculate frequency for each pathway
tfGO = zeros(length(path),length(pathGO));

for i=1:length(path)
    [junk ia ib] = intersect(path{i},pathGO);
    tfGO(i,ib) = 1;
end

pathGO(sum(tfGO)>1)'
inx = find(sum(tfGO)>1);

figure; 
imagesc(tfGO(:,inx)');
set(gca,'xtick',1:size(tfGO(:,inx),1),'xticklabel',screenLabels,'ytick',1:size(tfGO(:,inx),2),'yticklabel',pathGO(inx));

pathGO(inx)';

%% print the shared GO terms (used as input for REVIGO)
for i=1:length(inx)
    curLine = pathGO(inx(i));
    tokens = split(curLine);
    fprintf("%s\n",tokens{1});
end


%% calculate frequency for each pathway
tfKEGG = zeros(length(path),length(pathKEGG));

for i=1:length(path)
    [junk ia ib] = intersect(path{i},pathKEGG);
    tfKEGG(i,ib) = 1;
end

pathKEGG(sum(tfKEGG)>0)'
inx = sum(tfKEGG)>0;

figure; 
imagesc(tfKEGG(:,inx)');

set(gca,'xtick',1:size(tfKEGG(:,inx),1),'xticklabel',screenLabels,'ytick',1:size(tfKEGG(:,inx),2),'yticklabel',pathKEGG(inx));

%% plot by REVIGO
[num,txt,raw] = xlsread(inputFileName,'REVIGO (small)');
revigoCat = raw(:,1);
revigoDesc = raw(:,2);

figure; hold on;
for iCat = 1:length(revigoCat)
    curCat = revigoCat(iCat);
    for iScreen = 1:length(path)
        for iGOhit = 1:length(path{iScreen})
            curLine = path{iScreen}{iGOhit};
            tokens = split(curLine);
            curGOhit = tokens{1};
            if(strcmp(curCat{1},curGOhit))
                curQval = pval{iScreen}(iGOhit);
                curSetSize = setSize{iScreen}(iGOhit);
                curColor = min([-log10(curQval)/3 1]);
                plot(iScreen,iCat,'o','markersize',log10(curSetSize)*5,'markerfacecolor',[curColor 0 0],'markeredgecolor',[curColor 0 0]);
            end
        end
    end
end
set(gca,'xtick',1:length(screenLabels),'xticklabel',screenLabels,'ytick',1:length(revigoDesc),'yticklabel',revigoDesc);
set(gca,'xlim',[0.5 length(screenLabels)+0.5],'ylim',[0.5 length(revigoCat)+0.5]);
set(gcf,'position',[440   153   293   645]);
ytickangle(45); xtickangle(45);
box on; grid on;

% plot a "legend" image
fake_qval = [0.0005,0.005,0.05];
fake_setSize = [20,100,500];
figure; hold on;
for i=1:3
curQval = fake_qval(i);
                curSetSize = fake_setSize(i);
                curColor = min([-log10(curQval)/3 1]);
                plot(i,1,'o','markersize',log10(curSetSize)*5,'markerfacecolor',[curColor 0 0],'markeredgecolor',[curColor 0 0]);
end
set(gca,'xtick',1:length(screenLabels),'xticklabel',screenLabels,'ytick',1:length(revigoDesc),'yticklabel',revigoDesc);
set(gca,'xlim',[0.5 length(screenLabels)+0.5],'ylim',[0.5 length(revigoCat)+0.5]);
set(gcf,'position',[440   153   293   645]);
ytickangle(45); xtickangle(45);
box on; grid on;


    


