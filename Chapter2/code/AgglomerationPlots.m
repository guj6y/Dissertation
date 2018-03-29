clear
linkageType = 'Min';
%linkageType = 'Max';
load('data/Processed/webGeneration.mat')
T = load(sprintf('AgglomerationProps%sLinkage.mat',linkageType));

figureFont =  'CMU Serif';
propsLocal = T.propsLocal;
propsGlobal = T.propsGlobal;
localMeans = T.propsLocalMeans;
localCol = T.localCol;
globalCol = T.globalCol;
localMeanCol = T.localMeansCol;
fullProps = T.fullProps;
reducedProps = T.reducedProps;
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
webCodes = webCodes(webs);


jj = 0;

close all
localPropsPlotted = {
                        'vul'        ...
                        ,'gen'        ...
                        ,'patl'       ...
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

nLocalPropsPlotted = numel(localPropsPlotted);
globalPropsPlotted = {'C','Cf','Cff','Cfp','Cpf','Cpp','herb','int','omn','top','fPar','fCarn'};
    
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

selectThisLevelGlobal = arrayfun(@(x) find(propsGlobal(:,globalCol('web'))==x,1),webs);
Ss = propsGlobal(selectThisLevelGlobal,globalCol('S'));
Cs = propsGlobal(selectThisLevelGlobal,globalCol('C'));
              
alpha = 0.05;
m = nLocalPropsPlotted;
alphaSeq = alpha./(m+1-(1:m));

textColRejected = [0.8 0.2 0.2];
textColFailed = [0.5 0.5 0.5];
lineColRejected = [1 0.5 0.5];
lineColFailed = [0.8 0.8 0.8];

plotsToPlot = [false false false false true];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw plots of the initial webs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(1)
    
    localFigFirst = figure;
    
    
    Ss = propsGlobal(selectThisLevelGlobal,globalCol('S'));
    selectThisLevelLocal = arrayfun(@(x,y) (propsLocal(:,localCol('S'))==y)&(propsLocal(:,localCol('web'))==x),webs',Ss,'UniformOutput',false);
    initialWebs = sum(cell2mat(selectThisLevelLocal'),2)>0;
    initialWebsPropsLocal = propsLocal(initialWebs,:);

    initCarn = (initialWebsPropsLocal(:,localCol('free'))>0.5);
    initPara = (initialWebsPropsLocal(:,localCol('para'))>=0.5);
    initCarnOrPara = initCarn|initPara;

    cMapMin = round(min(propsGlobal(selectThisLevelGlobal,globalCol('C')))*.99-.001,3);
    cMapMax = round(max(propsGlobal(selectThisLevelGlobal,globalCol('C')))*1.01+.001,3);
    colMatrix = summer(100);
    cVals = linspace(cMapMin,cMapMax,100);

    
    tDiffs = zeros(nLocalPropsPlotted,1);
    pDiffs = zeros(nLocalPropsPlotted,1);
    xMeans = zeros(nLocalPropsPlotted,1);
    yMeans = zeros(nLocalPropsPlotted,1);
    
    for jj = 1:nLocalPropsPlotted
        figure(localFigFirst);
        xs = localMeans(selectThisLevelGlobal,localMeanCol(localPropsPlottedFree{jj}));
        ys = localMeans(selectThisLevelGlobal,localMeanCol(localPropsPlottedPara{jj}));
        gps = localMeans(selectThisLevelGlobal,localMeanCol('web'));
        markSizes = Ss/3;
        markColors = colMatrix(arrayfun(@(x)find(x<cVals,1),Cs),:);

        subplot(3,4,jj);
        h = gscatter(xs,ys,gps,markColors,'.',markSizes,'off');
        xlabel('Carnivore Average')
        ylabel('Parasite Average')
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
        pDiff = tcdf(abs(tDiff),numel(diffs)-1,'upper');
        pDiffs(jj) = pDiff;
        lb = min(xl(1),yl(1));
        ub = max(xl(2),yl(2));
        axis([lb ub lb ub]);
        rl = refline(1,0);
        rl.Color = [.9,.9,.9];
        uistack(h,'top')
    end
    
    [pDiffsSorted,Idx] = sort(pDiffs);
    [~,IdxIdx] = sort(Idx);
    rejected = false(nLocalPropsPlotted,1);
    rejected(Idx) = pDiffsSorted<alphaSeq';
    
    for jj = 1:nLocalPropsPlotted
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
    end
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
    set(localFigFirst, 'Units', 'Inches', 'Position', [0, 0, 17, 10])
    print(sprintf('../figures/initialProps%sLinkage.png',linkageType),'-dpng','-r0')
    
end

if plotsToPlot(2)||plotsToPlot(3)||plotsToPlot(4)
    allNodeLevels = propsGlobal(:,globalCol('S'));
    endings = [find(diff(allNodeLevels)>0); 
                numel(allNodeLevels)-1];
    beginnings = [1; endings(1:end-1)+1];
    endings = endings -5;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw plots of the agglomerated webs, size on x-axis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(2)
    localFigAllSizes = figure;
    
    xAllLevels = unique(allNodeLevels);
    nLevels = numel(xAllLevels);
    meanFAll = zeros(nLocalPropsPlotted,nLevels);
    meanPAll = zeros(nLocalPropsPlotted,nLevels);

    allTs = zeros(nLevels,nLocalPropsPlotted);
    allPs = zeros(nLevels,nLocalPropsPlotted);

    hEveryWeb = cell(nLocalPropsPlotted,nWebs);
    h = cell(nLocalPropsPlotted,1);
    hAll = cell(nLocalPropsPlotted,1);
    for jj = 1:nLocalPropsPlotted
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
    end
    
    rejectedAllLevels = zeros(nLevels,nLocalPropsPlotted);
    colsAllLevels = cell(nLevels,nLocalPropsPlotted);
    [colsAllLevels{:}] = deal(zeros(1,3));

    
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
    
    for jj = 1:nLocalPropsPlotted
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
        plot([Ss(end);Ss(end)],[yl(1) ;yl(2)],'Color',[.6,.6,.6])
        axis([0,max(Ss) ylim]);
        ylim(yl);
        hold off
        
    end
    
    arrayfun(@(x) set(x,'FontName',figureFont),localFigAllSizes.Children);
    set(localFigAllSizes, 'Units', 'Inches', 'Position', [0, 0, 17, 10])
    print(sprintf('../figures/allPropsBySize%sLinkage.png',linkageType),'-dpng','-r0')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This sets up some important distance stuff.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(3)||plotsToPlot(4)||plotsToPlot(5)
    allDistances = localMeans(:,localMeanCol('minDistance'));
    nEdges = 81;
    nBins = nEdges-1;
    dMin = 0;
    dMax = 0.4;
    distancesToPlot = linspace(dMin,dMax,nEdges);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw plots of the agglomerated, minimum node distance on x-axis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(3)
    localFigDistance = figure;
    hDistances = cell(nLocalPropsPlotted,nWebs);
    plotByMinDistance = nan(nEdges*nWebs,numel(localPropsPlotted));
    distanceGroups = repmat((1:nEdges)',nWebs,1);
    
    for ii = 1:nWebs
        distances_ii = allDistances(beginnings(ii):endings(ii));

        valuesFree_ii = localMeans(beginnings(ii):endings(ii),freeIndices);
        valuesPara_ii = localMeans(beginnings(ii):endings(ii),paraIndices); 
        diffs_ii = valuesFree_ii - valuesPara_ii;

        lastDistanceGreaterThan = sum(distances_ii<=distancesToPlot)';
        diffsToPlot_ii = diffs_ii(lastDistanceGreaterThan,:);

        plotByMinDistance(((1:nEdges)+nEdges*(ii-1)),:) = diffsToPlot_ii;

    end
    
    [groupedDistanceMeans, groupedDistanceStd, groupedDistanceN] = ...
        grpstats(plotByMinDistance, distanceGroups...
                    ,{   @(x) mean(x,'omitnan')...
                        ,@(x) std(x,'omitnan')...
                        ,@(x) sum(isfinite(x))...
                     }...
                );
            
    groupedDistanceTs = groupedDistanceMeans./(groupedDistanceStd./sqrt(groupedDistanceN));
    groupedDistancePs = tcdf(abs(groupedDistanceTs),groupedDistanceN-1,'upper');
    
    save(sprintf('../../Chapter3/code/groupedDistances%sLinkge.mat',linkageType),...
        'groupedDistanceMeans','groupedDistanceStd','groupedDistanceN'...
       ,'groupedDistanceTs','groupedDistancePs','plotByMinDistance'...
       ,'distancesToPlot');
   
   for jj = 1:nLocalPropsPlotted
   subplot(3,4,jj);
    hold on
        for ii = 1:nWebs
            thisWeb = (1:nEdges)+(ii-1)*nEdges;
            plotByMinDistanceThisWeb = plotByMinDistance(thisWeb,jj);
            hDistances{jj,ii} = plot(distancesToPlot,plotByMinDistanceThisWeb...
                ,'LineWidth',.5...
                ,'Color',[.7,.7,.7]);
        end
    hDistanceMean = plot(distancesToPlot,groupedDistanceMeans(:,jj)...
            ,'LineWidth',3 ...
            ,'Color',[0 0 1]);
    hold off
    title(fancyLocalNames{jj},'Interpreter','LaTeX')
    end
   
    rejectedAllBins = zeros(nEdges,nLocalPropsPlotted);
    colsAllBins = cell(nEdges,nLocalPropsPlotted);

    for ii = 1:nEdges
        PsThisBin = groupedDistancePs(ii,:);
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
   
    for jj = 1:nLocalPropsPlotted
        
        subplot(3,4,jj);
        hold on
        hDist = plot([distancesToPlot;distancesToPlot],[zeros(size(distancesToPlot)); groupedDistanceMeans(:,jj)']);
        set(hDist, {'color'}, colsAllBins(:,jj));
        set(hDist,'linewidth',3)
        uistack(hDistanceMean,'top')
        uistack(hDist,'bottom')
        rl = refline(0,0);
        rl.Color = [.5,.5,.5];
        rl.LineStyle = '--';
        yl = ylim();
        ylim(yl);
        xlabel('Min. Cluster Distance')
        ylabel('Carnivore $-$ Parasite','Interpreter','LaTeX')
        hold off

    end
    
    arrayfun(@(x) set(x,'FontName',figureFont),localFigDistance.Children);
    set(localFigDistance, 'Units', 'Inches', 'Position', [0, 0, 17, 10])
    print(sprintf('../figures/allPropsDist%sLinkage.png',linkageType),'-dpng','-r0')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw plots of the global Properties, minimum node distance on x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(4)
    globalPropsPlots = figure;
    for jj = 1:nLocalPropsPlotted
        subplot(3,4,jj);
        hold on
        for ii = 1:nWebs
            distances_ii = allDistances(beginnings(ii):endings(ii));

            globalPropsPlot = propsGlobal(beginnings(ii):endings(ii),globalCol(globalPropsPlotted{jj}));

            lastDistanceGreaterThan = sum(distances_ii<=distancesToPlot)';
            plotMe = globalPropsPlot(lastDistanceGreaterThan,:);


            hDistances{jj,ii} = plot(distancesToPlot,plotMe...
                ,'LineWidth',.5...
                ,'Color',[.7,.7,.7]);
        end
        hold off
        title(globalPropsPlotted{jj},'Interpreter','LaTeX')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw plots of the mean similarity between various classes of taxa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fullWebNames = {'(a) Bahia San Quintin'...
               ,'(b) Carpinteria Salt Marsh'...
               ,'(c) Estero de Punta Banda'...
               ,'(d) Otago Harbor'...
               ,'(e) Sylt Tidal Basin'...
               ,'(f) Flensburg Fjord'};
           
           shortWebNames = {'bahia','carp','punta','flens','otago','sylt','meta'};
if plotsToPlot(5)
    averageDistanceGroups = figure;

    useFull = true;
    for ii = 1:nWebs
        subplot(3,2,ii)
        title(webCodes{ii});
        if useFull
            usedProps = fullProps;
        else
            usedProps = reducedProps;
        end
                
        
        xPlot = cellfun(@(x) min(min(x)),usedProps{ii}.jac);
        paraParaJac = cellfun(@(s,t) triu(s(t>=0.5,t>=0.5),1),usedProps{ii}.jac(:),usedProps{ii}.para(:),'UniformOutput',false);
        paraCarnJac = cellfun(@(s,t,u) triu(s(t>=0.5,u>=0.5),1),usedProps{ii}.jac(:),usedProps{ii}.para(:),usedProps{ii}.carn(:),'UniformOutput',false);
        paraFreeJac = cellfun(@(s,t,u) triu(s(t>=0.5,u>=0.5),1),usedProps{ii}.jac(:),usedProps{ii}.para(:),usedProps{ii}.free(:),'UniformOutput',false);
        carnFreeJac = cellfun(@(s,t,u) triu(s(t>=0.5,u>=0.5),1),usedProps{ii}.jac(:),usedProps{ii}.carn(:),usedProps{ii}.free(:),'UniformOutput',false);
        freeFreeJac = cellfun(@(s,t) triu(s(t,t),1),usedProps{ii}.jac(:),usedProps{ii}.free(:),'UniformOutput',false);
        carnCarnJac = cellfun(@(s,t) triu(s(t>=0.5,t>=0.5),1),usedProps{ii}.jac(:),usedProps{ii}.carn(:),'UniformOutput',false);
        
            
        y1 = cellfun(@(x) mean(x(x>0)),paraParaJac);
        y2 = cellfun(@(x) mean(x(x>0)),paraCarnJac);
        y3 = cellfun(@(x) mean(x(x>0)),paraFreeJac);
        y4 = cellfun(@(x) mean(x(x>0)),freeFreeJac);
        y5 = cellfun(@(x) mean(x(x>0)),carnCarnJac);
        h = plot(xPlot,y1,'.',xPlot,y2,'.',xPlot,y3,'.',xPlot,y4,'.',xPlot,y5,'.');
        arrayfun(@(x) set(x,'MarkerSize',8),h);
        legend('Parasites','Par. to Carn.'...
            ,'Par. to F. Livers','Free Livers','Carnivores','Location','SouthEast')
        ylabel('Mean jaccard Distance')
        xlabel('Minimum Cluster Distance')
        title(fullWebNames{ii})
        arrayfun( @(x) set(x,'FontName','CMU Serif'),averageDistanceGroups.Children)
    end
    print(sprintf('../figures/meanClassDist.png'),'-dpng','-r0')
    
end
%print(sprintf('../figures/allPropsDist%sLinkage.png',linkageType),'-dpng','-r0')
