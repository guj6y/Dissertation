clear

load('data/Processed/webGeneration.mat')
T = load('AgglomerationProps.mat');

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
                        ...,'cc0'         ...
                        ,'pr'         ...
                        ,'btwns'       ...
                        ,'ecoBtwns'    ...
                        ,'meanVulPrey'...
                        ,'meanGenPred'...
                        ,'ccCyc'          ...
                        ,'ccMid'          ...
                        ,'ccIn'          ...
                        ,'ccOut'          ...
                        };
localPropsPlottedFree = cellfun(@(x) strcat(x,'Free'),localPropsPlotted,'UniformOutput',false);
localPropsPlottedPara = cellfun(@(x) strcat(x,'Para'),localPropsPlotted,'UniformOutput',false);

freeIndices = values(localMeanCol,localPropsPlottedFree);
freeIndices = [freeIndices{:}];

paraIndices = values(localMeanCol,localPropsPlottedPara);
paraIndices = [paraIndices{:}];

fancyLocalNames = {'(a) Vulnerability'...
                  ,'(b) Generality'...
                  ,'(c) Trophic Level'...
                  ,'(d) Eco. FlowRank'...
                  ,'(e) Betweenness Centrality'...
                  ,'(f) Eco. Betweenness Centrality'...
                  ,'(g) Mean Vul. Resources'...
                  ,'(h) Mean Gen. Consumers'...
                  ,'(i) Frac. 3-Cycles'...                  cycle
                  ,'(j) Frac. Resource-Consumer Links'...   mid
                  ,'(k) Frac. Resource Links'...            in
                  ,'(l) Frac. Consumer Links'};            %out

localFigFirst = figure;
localFigAll = figure;
localFigDistance = figure;
distanceSPlots = figure;
localFigRelNodeLevel = figure;
figureFont =  'CMU Serif';

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
allDistances = localMeans(:,localMeanCol('minDistance'));



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
hEveryWeb = cell(nLocalPropsPlotted,nWebs);
hDistances = cell(nLocalPropsPlotted,nWebs);

endings = [find(diff(allNodeLevels)>0);numel(allNodeLevels)-1];
beginnings = [1; endings(1:end-1)+1];
endings = endings -5;
nEdges = 21;
nBins = nEdges-1;
dMin = 0;
dMax = 0.8;
distancesToPlotEdges = linspace(dMin,dMax,nEdges);
binnedValues = nan(nBins*nWebs,numel(localPropsPlotted));
deltaDistance = mean(diff(distancesToPlotEdges));
distancesToPlot = linspace(dMin+deltaDistance/2,dMax-deltaDistance/2,nBins);
grpsForDistancePlot = repmat((1:nBins)',nWebs,1);
for ii = 1:nWebs
    distances_ii = allDistances(beginnings(ii):endings(ii));
    
    valuesFree_ii = localMeans(beginnings(ii):endings(ii),freeIndices);
    valuesPara_ii = localMeans(beginnings(ii):endings(ii),paraIndices); 
    diffs_ii = valuesFree_ii - valuesPara_ii;
    discretizedDistances_ii = discretize(distances_ii,distancesToPlotEdges);
    distanceDiscretizedMeans_ii = grpstats(diffs_ii,discretizedDistances_ii);
    discretizedDistances_ii(isnan(discretizedDistances_ii)) = [];
    discretizedDistances_ii = unique(discretizedDistances_ii);
    binnedValues(discretizedDistances_ii+(nBins*(ii-1)),:) = distanceDiscretizedMeans_ii;
    
end
nanBinnedValues = isnan(binnedValues(:,1));
[meansToPlotBinnedDistances,stdsToPlotBinnedDistances,nsToPlotBinnedDistances]...
    = grpstats(binnedValues,grpsForDistancePlot...
        ,{@(x) mean(x,'omitnan'),@(x) std(x,'omitnan'),@numel});
TsBinnedDistances = meansToPlotBinnedDistances./(stdsToPlotBinnedDistances./sqrt(nsToPlotBinnedDistances));
PsBinnedDistances = tcdf(abs(TsBinnedDistances),nsToPlotBinnedDistances-1,'upper');

%This code plots distance versus fraction of S remaining. shows very high
%correlation between the two.
    %{
figure(distanceSPlots);
    for ii = 1:nWebs
        hold on
        
        plot(localMeans(beginnings(ii):endings(ii),localMeanCol('minDistance')),allNodeLevels(beginnings(ii):endings(ii))/allNodeLevels(beginnings(ii))...
            ,'Color',[.5,.5,.5]...
            );
        hold off
    end
    q = refline(-1,1);
    q.LineStyle = '--';
%}
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
    ts = meanD_jj./(stdD_jj./sqrt(numelD_jj));
    ts(numelD_jj==1) = nan;
    Ps = tcdf(abs(ts),numelD_jj-1,'upper');
    
    allTs(:,jj) = ts;
    allPs(:,jj) = Ps;
    subplot(3,4,jj);
    
    hAll{jj} = plot(xAllLevels,meanF_jj,'b',xAllLevels,meanP_jj,'m');
    set(hAll{jj},{'LineWidth'},{1;1});
    
    for ii = 1:nWebs
        hold on
        hEveryWeb{jj,ii} = plot(allNodeLevels(beginnings(ii):endings(ii)),freesThisProperty(beginnings(ii):endings(ii))...
            ,allNodeLevels(beginnings(ii):endings(ii)),parasThisProperty(beginnings(ii):endings(ii))...
            );
        set(hEveryWeb{jj,ii},{'LineWidth'},{.5;.5})
        set(hEveryWeb{jj,ii},{'Color'},{[.7,.7,1];[1 .7 1]})
        hold off
        
    end
    
    title(fancyLocalNames{jj},'Interpreter','LaTeX')
    xlabel('S')
    ylabel('Mean')
    
    figure(localFigDistance) 
    subplot(3,4,jj);
    hold on
        for ii = 1:nWebs
            thisWeb = (1:nBins)+(ii-1)*nBins;
            binnedValuesThisWeb = binnedValues(thisWeb,jj);
            nanBinnedValuesThisWeb = nanBinnedValues(thisWeb);
            grpsForDistancePlotThisWeb = nanBinnedValues(thisWeb);
            hDistances{jj,ii} = plot(distancesToPlot(~nanBinnedValuesThisWeb),binnedValuesThisWeb(~nanBinnedValuesThisWeb)...
                ,'LineWidth',.5...
                ,'Color',[.7,.7,.7]);
        end
    hDistanceMean = plot(distancesToPlot,meansToPlotBinnedDistances(:,jj)...
            ,'LineWidth',3 ...
            ,'Color',[0 0 1]);
    hold off
    title(fancyLocalNames{jj},'Interpreter','LaTeX')
    
    figure(localFigRelNodeLevel) 
    for ii = 1:nWebs
        subplot(3,4,jj);
        hold on
        %{
        hEveryWeb{jj,ii} = plot(localMeans(beginnings(ii):endings(ii),localMeanCol('minDistance')),freesThisProperty(beginnings(ii):endings(ii))...
            ,localMeans(beginnings(ii):endings(ii),localMeanCol('minDistance')),parasThisProperty(beginnings(ii):endings(ii))...
            );
        set(hEveryWeb{jj,ii},{'LineWidth'},{.5;.5})
        set(hEveryWeb{jj,ii},{'Color'},{[.7,.7,1];[1 .7 1]})
        %}
        q = plot(1-allNodeLevels(beginnings(ii):endings(ii))/allNodeLevels(beginnings(ii)),freesThisProperty(beginnings(ii):endings(ii))...
            - parasThisProperty(beginnings(ii):endings(ii))...
            );
        if ii == 4
            q.Color = 'r';
        else
            q.Color =[.5,.5,.5];
        end
        hold off
    end
    q = refline(0,0);
    q.LineStyle = '--';
    q.Color = 'k';
    title(fancyLocalNames{jj},'Interpreter','LaTeX')
    
    figure(localFigFirst);
    xs = localMeans(selectThisLevelGlobal,localMeanCol(localPropsPlottedFree{jj}));
    ys = localMeans(selectThisLevelGlobal,localMeanCol(localPropsPlottedPara{jj}));
    gps = localMeans(selectThisLevelGlobal,localMeanCol('web'));
    markSizes = Ss/3;
    markColors = colMatrix(arrayfun(@(x)find(x<cVals,1),Cs),:);
    
    subplot(3,4,jj);
    h = gscatter(xs,ys,gps,markColors,'.',markSizes,'off');
    title(fancyLocalNames{jj},'Interpreter','LaTeX')
    if jj == 8
        xlabel('Carnivore Average')
    end
    if jj == 4
    
        ylabel('Parasitic Average')
    end
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
m = nLocalPropsPlotted;
alphaSeq = alpha./(m+1-(1:m));
[pDiffsSorted,Idx] = sort(pDiffs);
[~,IdxIdx] = sort(Idx);
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

rejectedAllBins = zeros(nBins,nLocalPropsPlotted);
colsAllBins = cell(nBins,nLocalPropsPlotted);

for ii = 1:nBins
    PsThisBin = PsBinnedDistances(ii,:);
    [pSorted,Idx] = sort(PsThisBin);
    thisBinRejected = false(nLocalPropsPlotted,1);
    thisBinRejected(Idx) = pSorted<alphaSeq;
    if sum(isnan(pSorted))>0
        [colsAllBins{ii,:}] = deal('none');
    else
        [colsAllBins{ii,thisBinRejected}] = deal(lineColRejected);
        [colsAllBins{ii,~thisBinRejected}] = deal(lineColFailed);
    end
    
    
end
textColRejected = [0.8 0.2 0.2];
textColFailed = [0.5 0.5 0.5];
for jj = 1:nLocalPropsPlotted
        figure(localFigFirst);
        subplot(3,4,jj);
        hold on
        xl = xlim();
        yl = ylim();
        
        if rejected(jj)
            textCol = textColRejected;
            plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
            plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
            meanLabel = sprintf('$t=%.3g$\n$P=%.3g<%.3g$',tDiffs(jj),pDiffs(jj),alphaSeq(IdxIdx(jj)));
        else
            textCol = textColFailed;
            plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
            plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
            meanLabel = sprintf('$t=%.3g$\n$P=%.3g\\geq%.3g$',tDiffs(jj),pDiffs(jj),alphaSeq(IdxIdx(jj)));
        end
        
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
                        ,'Interpreter','LaTeX'...
                         );
        hold off
        axis([xl yl]);
        
        figure(localFigDistance)
        subplot(3,4,jj);
        hold on
        h = plot([distancesToPlot;distancesToPlot],[zeros(size(distancesToPlot)); meansToPlotBinnedDistances(:,jj)']);
        set(h, {'color'}, colsAllBins(:,jj));
        set(h,'linewidth',3)
        uistack(hDistanceMean,'top')
        uistack(h,'bottom')
        rl = refline(0,0);
        rl.Color = [.5,.5,.5];
        rl.LineStyle = '--';
        yl = ylim();
        ylim(yl);
        hold off
        
        figure(localFigAll);
        subplot(3,4,jj);
        hold on
        h = plot([xAllLevels';xAllLevels'],[meanFAll(jj,:);meanPAll(jj,:)]);
        set(h, {'color'}, colsAllLevels(:,jj));
        set(h,'linewidth',3)
        uistack(h,'bottom')
        cellfun(@(x)uistack(x,'top'),hEveryWeb(jj,:),'UniformOutput',false);
        uistack(hAll{jj},'top')
        rl = refline(0,1);
        rl.Color = [.5,.5,.5];
        rl.LineStyle = '--';
        yl = ylim();
        plot([Ss(end);Ss(end)],[yl(1);yl(2)],'Color',[.6,.6,.6])
        ylim(yl);
        hold off
end
    %}
%

figure(localFigFirst);
hp9 = get(subplot(3,4,12),'Position');
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
set(localFigFirst, 'Units', 'Inches', 'Position', [0, 0, 17, 10])
set(localFigAll, 'Units', 'Inches', 'Position', [0, 0, 17, 10])
set(localFigDistance, 'Units', 'Inches', 'Position', [0, 0, 17, 10])
%
figure(localFigFirst);
print('../figures/initialProps.png','-dpng','-r0')
figure(localFigAll);
print('../figures/allProps.png','-dpng','-r0')
figure(localFigDistance);
print('../figures/allPropsDist.png','-dpng','-r0')
%}
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
