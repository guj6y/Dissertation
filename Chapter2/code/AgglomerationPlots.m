clear

load('data/Processed/webGeneration.mat')
T = load('AgglomerationProps.mat');

figureFont =  'CMU Serif';
propsLocal = T.propsLocal;
propsGlobal = T.propsGlobal;
localMeans = T.propsLocalMeans;
localCol = T.localCol;
globalCol = T.globalCol;
localMeanCol = T.localMeansCol;
clear T

[nSpecies,nLocalVars] = size(propsLocal);
[nClusters,nGlobalVars] = size(propsGlobal);
[~,nLocalMeans] = size(localMeans);

webs = 1:6; %6 empirical webs + averages for nichewebs with same S,C as these webs.
nWebs = numel(webs);


%globalFig0 = figure();
%localFig0 = figure();
xmin = inf;
ymin = inf;
xmax = 0;
ymax = 0;
webs = [3 2 6 1 5 4];

allYs = zeros(nWebs,9);
allXs = zeros(nWebs,9);

webCodes = {'BSQ','CSM','EPB','FF','OH','STB'};
varNames = {'$VUL$','$GEN$','$PATL$','$CC$','$EIG$','$C_B$','$C_EB$','$\overline{v_r}$','$\overline{g_c}$'};
webCodes = webCodes(webs);


jj = 0;

close all
localPropsPlotted = {
                        'vul'        ...
                        ,'gen'        ...
                        ,'patl'       ...
                        ,'cc'         ...
                        ,'pr'         ...
                        ,'btwns'       ...
                        ,'ecoBtwns'    ...
                        ,'meanVulPrey'...
                        ,'meanGenPred'...
                        };
localPropsPlottedFree = cellfun(@(x) strcat(x,'Free'),localPropsPlotted,'UniformOutput',false);
localPropsPlottedPara = cellfun(@(x) strcat(x,'Para'),localPropsPlotted,'UniformOutput',false);

fancyLocalNames = {'(a) v','(b) g','(c) TL','(d) cc','(e) \lambda','(f) C_B','(g) C_{EB}','(h) \langle v_r\rangle','(i) \langle g_c\rangle'};

localFigFirst = figure;
localFigAll = figure;

webCount = 0;
plottingInitial = true;
initialSizes = propsGlobal(arrayfun(@(x) find(propsGlobal(:,globalCol('web'))==x,1),webs),globalCol('S'));
if plottingInitial
    selectThisLevelGlobal = arrayfun(@(x) find(propsGlobal(:,globalCol('web'))==x,1),webs);
else
    plotAtThisLevel = 50;
    selectThisLevelLocal = propsLocal(:,localCol('S')) == plotAtThisLevel;
    selectThisLevelGlobal = propsGlobal(:,globalCol('S')) == plotAtThisLevel;
end

Ss = propsGlobal(selectThisLevelGlobal,globalCol('S'));
Cs = propsGlobal(selectThisLevelGlobal,globalCol('C'));
cMapMin = round(min(propsGlobal(selectThisLevelGlobal,globalCol('C')))*.99-.001,3);
cMapMax = round(max(propsGlobal(selectThisLevelGlobal,globalCol('C')))*1.01+.001,3);
colMatrix = summer(100);
cVals = linspace(cMapMin,cMapMax,100);

allNodeLevels = propsGlobal(:,globalCol('S'));
xAllLevels = unique(allNodeLevels);
nLocalPropsPlotted = numel(localPropsPlotted);
nLevels = numel(xAllLevels);

meanFAll = zeros(nLocalPropsPlotted,nLevels);
meanPAll = zeros(nLocalPropsPlotted,nLevels);

tDiffs = zeros(nLocalPropsPlotted,1);
pDiffs = zeros(nLocalPropsPlotted,1);
xMeans = zeros(nLocalPropsPlotted,1);
yMeans = zeros(nLocalPropsPlotted,1);

allTs = zeros(nLevels,nLocalPropsPlotted);
allPs = zeros(nLevels,nLocalPropsPlotted);

h = cell(nLocalPropsPlotted,1);
hAll = cell(nLocalPropsPlotted,1);
for jj = 1:nLocalPropsPlotted
    figure(localFigAll);
    freesThisProperty = localMeans(:,localMeanCol(localPropsPlottedFree{jj}));
    parasThisProperty = localMeans(:,localMeanCol(localPropsPlottedPara{jj}));
    diffsThisProperty = freesThisProperty-parasThisProperty;
    [meanF_jj] = grpstats(freesThisProperty,allNodeLevels);
    [meanP_jj] = grpstats(parasThisProperty,allNodeLevels);
    meanFAll(jj,:) = meanF_jj;
    meanPAll(jj,:) = meanP_jj;
    [meanD_jj,numelD_jj,stdD_jj] = grpstats(diffsThisProperty,allNodeLevels,{@mean,@numel,@std});
    ts = meanD_jj./stdD_jj;
    ts(numelD_jj==1) = nan;
    Ps = tcdf(abs(ts),numelD_jj-1,'upper');
    
    allTs(:,jj) = ts;
    allPs(:,jj) = Ps;
    subplot(3,3,jj);
    
    hAll{jj} = plot(xAllLevels,meanF_jj,'b',xAllLevels,meanP_jj,'m');
    set(hAll{jj},{'LineWidth'},{1;1});
    title(fancyLocalNames{jj})
    xlabel('S')
    ylabel('Mean')
    
    
    figure(localFigFirst);
    xs = localMeans(selectThisLevelGlobal,localMeanCol(localPropsPlottedFree{jj}));
    ys = localMeans(selectThisLevelGlobal,localMeanCol(localPropsPlottedPara{jj}));
    gps = localMeans(selectThisLevelGlobal,localMeanCol('web'));
    markSizes = Ss/3;
    markColors = colMatrix(arrayfun(@(x)find(x<cVals,1),Cs),:);
    
    subplot(3,3,jj);
    h = gscatter(xs,ys,gps,markColors,'.',markSizes,'off');
    title(fancyLocalNames{jj},'Interpreter','tex')
    xlabel('Free-Livers')
    ylabel('Parasites')
    xl = xlim;
    yl = ylim;
    xMeans(jj) = mean(xs);
    yMeans(jj) = mean(ys);
    diffs = ys-xs;
    meanDiff = mean(diffs);
    stdDiff = std(diffs)/sqrt(numel(diffs));
    tDiff = meanDiff./stdDiff;
    tDiffs(jj) = tDiff;
    pDiff = tcdf(abs(tDiff),numel(tDiffs)-1,'upper');
    pDiffs(jj) = pDiff;
    lb = min(xl(1),yl(1));
    ub = max(xl(2),yl(2));
    axis([lb ub lb ub]);
    rl = refline(1,0);
    rl.Color = [.9,.9,.9];
    uistack(h,'top')
        
end
    
    
alpha = 0.05;
m = 9;
alphaSeq = alpha./(m+1-(1:m));
[pDiffsSorted,Idx] = sort(pDiffs);
rejected = false(nLocalPropsPlotted,1);
rejected(Idx) = pDiffsSorted<alphaSeq';
%Replicate above for each level.
rejectedAllLevels = zeros(nLevels,nLocalPropsPlotted);
colsAllLevels = cell(nLevels,nLocalPropsPlotted);
[colsAllLevels{:}] = deal(zeros(1,3));
lineColRejected = [1 0.5 0.5];
lineColFailed = [0.8 0.8 0.8];
for ii = 1:nLevels
    
    PsThisLevel = allPs(ii,:);
    [pSorted,Idx] = sort(PsThisLevel);
    rejectedThisLevel = false(nLocalPropsPlotted,1);
    rejectedThisLevel(Idx) = pSorted<alphaSeq;
    if sum(isnan(pSorted))>0
        [colsAllLevels{ii,:}] = deal('none');
    else
        [colsAllLevels{ii,rejectedThisLevel}] = deal(lineColRejected);
        [colsAllLevels{ii,~rejectedThisLevel}] = deal(lineColFailed);
    end
    
end
textColRejected = [0.8 0.2 0.2];
textColFailed = [0.5 0.5 0.5];
for jj = 1:9
        figure(localFigFirst);
        subplot(3,3,jj);
        hold on
        xl = xlim();
        yl = ylim();
        
        if rejected(jj)
            textCol = textColRejected;
            plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
            plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
        else
            textCol = textColFailed;
            plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
            plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
        end
        meanLabel = sprintf('t=%.3g\nP=%.3g',tDiffs(jj),pDiffs(jj));
        textLocOpts = [xl(1)+diff(xl)*.02,yl(1)+diff(yl)*.02;...LowerLeft
                       xl(1)+diff(xl)*.02,yl(2)-diff(yl)*.02;...UpperLeft
                       xl(2)-diff(xl)*.02,yl(2)-diff(yl)*.02;...UpperRight
                       xl(2)-diff(xl)*.02,yl(1)+diff(yl)*.02]; %LowerRight
                   
        alignmentOpts = {'VerticalAlignment','bottom','HorizontalAlignment','left';
                        'VerticalAlignment','top','HorizontalAlignment','left';
                        'VerticalAlignment','top','HorizontalAlignment','right';
                        'VerticalAlignment','bottom','HorizontalAlignment','right'};
        distToCentroids = sqrt(sum((textLocOpts-[xMeans(jj),yMeans(jj)]).^2,2));
        thisLoc = distToCentroids==max(distToCentroids);
        textLoc = textLocOpts(thisLoc,:);
        thisAlign = alignmentOpts(thisLoc,:);
        text(textLoc(1),textLoc(2),meanLabel...
                        ,'Color',textCol...
                        ,'FontSize',8 ...
                        ,thisAlign{:}...
                        ,'Margin',1 ...
                        ,'backgroundColor','white'...
                        ,'FontName',figureFont...
                         );
        hold off
        axis([xl yl]);
        
        figure(localFigAll);
        subplot(3,3,jj);
        hold on
        h = plot([xAllLevels';xAllLevels'],[meanFAll(jj,:);meanPAll(jj,:)]);
        set(h, {'color'}, colsAllLevels(:,jj));
        set(h,'linewidth',2)
        uistack(hAll{jj},'top')
        rl = refline(0,1);
        rl.Color = [.5,.5,.5];
        rl.LineStyle = '--';
        yl = ylim();
        plot([Ss(end);Ss(end)],[yl(1);yl(2)],'Color',[.1,.1,.1])
        ylim(yl);
        hold off
end
    %}
%

figure(localFigFirst);
hp9 = get(subplot(3,3,9),'Position');
a = colorbar('Position', [0.92 0.11 0.02 .81]...
             ,'Limits',[cMapMin,cMapMax]...
             );

cMapMid = round(mean([cMapMin,cMapMax]),3);
if min(abs(cMapMid-Cs))<.002
    cMapMid = [];
end
newTicks = [cMapMin; cMapMax; cMapMid; Cs];
[newTicks, Idx] = sort(newTicks);
a.Ticks = newTicks;
webCodes = webCodes';
if isempty(cMapMid)
    newTickLabels = [{cMapMin;cMapMax};webCodes(:)];  
else
    
    newTickLabels = [{cMapMin;cMapMax;cMapMid};webCodes(:)];
end

a.TickLabels = newTickLabels(Idx);
caxis([cMapMin,cMapMax])
title(a,'C')
colormap(summer);
arrayfun(@(x) set(x,'FontName',figureFont),localFigFirst.Children);
arrayfun(@(x) set(x,'FontName',figureFont),localFigAll.Children);
set(localFigFirst, 'Units', 'Inches', 'Position', [0, 0, 6, 6])
set(localFigAll, 'Units', 'Inches', 'Position', [0, 0, 6, 6])

figure(localFigFirst);
print('../figures/initialProps.png','-dpng','-r0')
figure(localFigAll);
print('../figures/allProps.png','-dpng','-r0')

%}
%{
for ii = 1:nWebs
    globalFigAll = figure();
    x = props{ii}.(linkageType).mean.global.nNodes; 
propCount = 0;    
 for jj = [3,8,4,11,5,12,6,13,7,14]
        
        propCount = propCount+1;
        subplot(5,2,propCount)
        hold on
            title(globalNames{jj})
        plot(x,props{ii}.(linkageType).mean.global.(globalNames{jj}),'b')
        plot(x,props{ii+6}.(linkageType).mean.global.(globalNames{jj}),'k')
        hold off
 end
end
%}
