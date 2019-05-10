%clear
linkageType = 'Max';
%linkageType = 'Max';
load -v7 'data/Processed/Generation.mat'
T = load(sprintf('AgglomerationProps%sLinkage.mat',linkageType));

localNames = {'web'           ...1
                ,'vul'        ...2
                ,'gen'        ...3
                ,'patl'       ...4
                ,'cc'         ...5
                ,'pr'         ...6
                ,'btwn'       ...7
                ,'ecoBtwn'    ...8
                ,'meanVulPrey'...9
                ,'meanGenPred'...10
                ,'para'       ...11
                ,'free'       ...12
                ,'S'          ...13
                ,'C'          ...14
                ,'size'       ...15
                ,'ccCyc'      ...16
                ,'ccMid'      ...17
                ,'ccIn'       ...18
                ,'ccOut'      ...19
                ,'counter'    ...20
                };


globalNames = { 'web'           ...1
              ,'S'              ...2
              ,'Cf'             ...3
              ,'Cff'            ...4
              ,'Cpp'            ...5
              ,'Cfp'            ...6
              ,'Cpf'            ...7
              ,'fPar'           ...8
              ,'numLargest'     ...9
              ,'n'              ...10
              ,'top'            ...11
              ,'int'            ...12
              ,'bas'            ...13
              ,'herb'           ...14
              ,'omn'            ...15
              ,'C'              ...16
              ,'d_J'            ...17
              ,'fCarn'          ...18
              ,'Cp'             ...19
              ,'Sf'             ...20
              ,'counter'        ...
              };
localMeanNames = {
                 'web'...                     1
                ,'vulFree'...                 2
                ,'vulPara'...                 3
                ,'genFree'...                 4
                ,'genPara'...                 5
                ,'patlFree'...                6
                ,'patlPara'...                7
                ,'ccFree'...                  8
                ,'ccPara'...                  9
                ,'prFree'...                  10
                ,'prPara'...                  11
                ,'btwnFree'...               12
                ,'btwnPara'...               13
                ,'ecoBtwnFree'...            14
                ,'ecoBtwnPara'...            15
                ,'meanVulPreyFree'...         16
                ,'meanVulPreyPara'...         17
                ,'meanGenPredFree'...         18
                ,'meanGenPredPara'...         19
                ,'S'...                       20
                ,'C'...                       21
                ,'minDistance'...             22
                ,'ccCycFree'        ...       23
                ,'ccCycPara'        ...       24
                ,'ccMidFree'        ...       25
                ,'ccMidPara'        ...       26
                ,'ccInFree'        ...        27
                ,'ccInPara'        ...        28
                ,'ccOutFree'        ...       29
                ,'ccOutPara'        ...       30
                ,'counter'                 ...31
            };
localCol = containers.Map(localNames, 1:numel(localNames));
globalCol = containers.Map(globalNames, 1:numel(globalNames));
localMeansCol = containers.Map(localMeanNames, 1:numel(localMeanNames));
    
figureFont =  'Times New Roman';
propsLocal = T.propsLocal;
propsGlobal = T.propsGlobal;
localMeans = T.propsLocalMeans;
fullProps = T.fullProps;
reducedProps = T.reducedProps;
clear T

[nSpecies,nLocalVars] = size(propsLocal);
[nClusters,nGlobalVars] = size(propsGlobal);
[~,nLocalMeans] = size(localMeans);

webs = 1:6; %6 empirical webs + averages for nichewebs with same S,C as these webs.
nWebs = numel(webs);


%globalFig0 = figure(); %localFig0 = figure();
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
localPropsPlotted = {'vul'         ...
                     ,'gen'        ...
                     ,'patl'       ...
                     ,'pr'         ...
                     ,'btwn'      ...
                     ,'ecoBtwn'   ...
                     ,'meanVulPrey'...
                     ,'meanGenPred'...
                     ,'ccOut'      ...consumer
                     ,'ccIn'       ...resource
                     ,'ccMid'      ...con-res
                     ,'ccCyc'      ...res-con
                     };

fancyLocalNames = {'(a) Vulnerability'...
                  ,'(b) Generality'...
                  ,'(c) Trophic Level'...
                  ,'(d) Eco. FlowRank'...
                  ,'(e) Betweenness Centrality'...
                  ,'(f) Eco. Betweenness Centrality'...
                  ,'(g) Mean Vul. Resources'...
                  ,'(h) Mean Gen. Consumers'...
                  ,'(i) Consumer Clustering'...CC_Out
                  ,'(j) Resource Clustering'...CC_In
                  ,'(k) Con. - Res. Clustering'...CC_Mid
                  ,'(l) Res. - Con. Clustering'...CC_Cyc
                  }; 
                        
nLocalPropsPlotted = numel(localPropsPlotted);
globalPropsPlotted = {'C','Cf','Cff','Cfp','Cpf','Cpp','herb','int','omn','top','fPar','fCarn'};
    
localPropsPlottedFree = cellfun(@(x) strcat(x,'Free'),localPropsPlotted,'UniformOutput',false);
localPropsPlottedPara = cellfun(@(x) strcat(x,'Para'),localPropsPlotted,'UniformOutput',false);

% freeIndices = values(localMeansCol,localPropsPlottedFree);
% freeIndices = [freeIndices{:}];
freeIndices = cellfun(@(x) localMeansCol(x), localPropsPlottedFree);

% paraIndices = values(localMeansCol,localPropsPlottedPara);
% paraIndices = [paraIndices{:}];
paraIndices = cellfun(@(x) localMeansCol(x), localPropsPlottedPara);

          

selectThisLevelGlobal = arrayfun(@(x) find(propsGlobal(:,globalCol('web'))==x,1),webs);
Ss = propsGlobal(selectThisLevelGlobal,globalCol('S'));
Cs = propsGlobal(selectThisLevelGlobal,globalCol('C'));
              
alpha = 0.05;
m = nLocalPropsPlotted;
alphaSeq = alpha./(m+1-(1:m));

textColRejected = [0.9 0.2 0.2];
textColFailed = [0.1 0.1 0.3];
lineColRejected = [1 0.5 0.5];
lineColFailed = [0.8 0.8 0.8];

plotsToPlot = [false true false false false false false false];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw plots of the initial webs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(1)
    
    AR = 9/5.5;
    img_width = 18;
    img_height = img_width/AR;
    # Setup geometry of figure:
    left_offset = 0.5;
    bot_offset = 0.04;
    h_fig_sep = .02;
    v_fig_sep = 0.02;
    v_rem = 1 - bot_offset - 3*v_fig_sep;
    h_rem = 1 - .11 - h_fig_sep*3 - left_offset;
    fig_width = h_rem/4;
    fig_heighth = v_rem/3;
    
    localFigFirst = figure('units','inches'...
                           ,'position',[10,10,img_width,img_height]); 
                           
    colormap(summer);
    cm_ = get(localFigFirst,'colormap');
    set(localFigFirst,'colormap',cm_(size(cm_,1):-1:1,:));
    
    
    Ss = propsGlobal(selectThisLevelGlobal,globalCol('S'));
    selectThisLevelLocal = arrayfun(@(x,y) (propsLocal(:,localCol('S'))==y)&(propsLocal(:,localCol('web'))==x),webs',Ss,'UniformOutput',false);
    initialWebs = sum(cell2mat(selectThisLevelLocal'),2)>0;
    initialWebsPropsLocal = propsLocal(initialWebs,:);

    initCarn = (initialWebsPropsLocal(:,localCol('free'))>0.5);
    initPara = (initialWebsPropsLocal(:,localCol('para'))>=0.5);
    initCarnOrPara = initCarn|initPara;

    cMapMin = min(propsGlobal(selectThisLevelGlobal,globalCol('C')))*.99-.001;
    cMapMax = max(propsGlobal(selectThisLevelGlobal,globalCol('C')))*1.01+.001;
    colMatrix = summer(100);    
    colMatrix(size(colMatrix,1):-1:1,:);
    cVals = linspace(cMapMin,cMapMax,100);

    
    tDiffs = zeros(nLocalPropsPlotted,1);
    pDiffs = zeros(nLocalPropsPlotted,1);
    xMeans = zeros(nLocalPropsPlotted,1);
    yMeans = zeros(nLocalPropsPlotted,1);
    
    rC = zeros(nLocalPropsPlotted,1);
    rS = zeros(nLocalPropsPlotted,1);
    
    for jj = 1:nLocalPropsPlotted
        figure(localFigFirst);
        xs = localMeans(selectThisLevelGlobal,localMeansCol(localPropsPlottedFree{jj}));
        ys = localMeans(selectThisLevelGlobal,localMeansCol(localPropsPlottedPara{jj}));
        gps = localMeans(selectThisLevelGlobal,localMeansCol('web'));
        markSizes = sqrt(Ss-80)*30;
        markColors = colMatrix(arrayfun(@(x)find(x<cVals,1),Cs),:);
        fig_col = mod(jj-1,4) + 1;
        fig_row = floor((jj-1)/4)+1;
        h_pos = left_offset + (fig_row-1)*(fig_width + h_fig_sep)
        v_pos = bot_offset + (fig_col-1)*(fig_heighth + v_fig_sep)
        subplot(3,4,jj,'position',[h_pos, v_pos, fig_width, fig_heighth]);
        h = scatter(xs,ys,markSizes,markColors,'filled');
        rC(jj) = corr(ys-xs,Cs);
        rS(jj) = corr(ys-xs,Ss);
        if jj > 8
            xlabel('Carnivore Average')
        else
            xlabel('')
        end
        if mod(jj,4) == 1
            ylabel('Parasite Average')
        else
            ylabel('')
        end
        title(fancyLocalNames{jj}...
            ,'Interpreter','tex'...
            ,'FontSize',14)
        xl = xlim;
        yl = ylim;
        xMeans(jj) = mean(xs);
        yMeans(jj) = mean(ys);
        diffs = ys-xs;
        meanDiff = mean(diffs);
        stdDiff = std(diffs)/sqrt(numel(diffs));
        tDiff = meanDiff./stdDiff;
        tDiffs(jj) = tDiff;
        pDiff = 1-tcdf(abs(tDiff),numel(diffs)-1);
        pDiffs(jj) = pDiff;
        lb = min(xl(1),yl(1));
        ub = max(xl(2),yl(2));
        
        %rl = refline(1,0);
        %rl.Color = textColFailed;
        % uistack(h,'top')
        hold on

        rl = line([lb ub],[lb ub],'color',textColFailed);
        h = scatter(xs,ys,markSizes,markColors);
        axis([lb ub lb ub]);
    end
    
    [pDiffsSorted,Idx] = sort(pDiffs);
    [~,IdxIdx] = sort(Idx);
    rejected = false(nLocalPropsPlotted,1);
    rejected(Idx) = pDiffsSorted<alphaSeq';
    
    tC = rC*2./sqrt(1-rC.^2);
    tS = rS*2./sqrt(1-rS.^2);
    pC = (1-tcdf(abs(tC),4))*2;
    pS = (1-tcdf(abs(tS),4))*2;
    
    CS_tab = [rC tC pC rS tS pS];
    data_cell = horzcat(fancyLocalNames',mat2cell(CS_tab,repmat([1],12,1),repmat([1],6,1)));
    fid = fopen('interaction_cs.csv','w')
    fprintf(fid,',corr(C),t(C),P(C),corr(S),t(S),P(S)\n')
    fprintf(fid,'%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n',data_cell'{:})
    fclose(fid)
    for jj = 1:nLocalPropsPlotted
       fig_col = mod(jj-1,4) + 1;
       fig_row = floor((jj-1)/4)+1;
       h_pos = left_offset + (fig_row-1)*(fig_width + h_fig_sep);
       v_pos = bot_offset + (fig_col-1)*(fig_heighth + v_fig_sep);
       subplot(3,4,jj,'position',[h_pos, v_pos, fig_width, fig_heighth]);
        hold on
        xl = xlim();
        yl = ylim();
        
        if rejected(jj)
            textCol = textColRejected;
            plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
            plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
            meanLabel = sprintf('t=%.3g\nP=%.3g<%.3g',tDiffs(jj),pDiffs(jj),alphaSeq(IdxIdx(jj)));
        else
            textCol = textColFailed;
            plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
            plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
            meanLabel = sprintf('t=%.3g\nP=%.3g\\geq%.3g',tDiffs(jj),pDiffs(jj),alphaSeq(IdxIdx(jj)));
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
                        ,'FontSize',14 ...
                        ,thisAlign{:}...
                        ,'Margin',1 ...
                        ,'backgroundColor','white'...
                        ,'FontName',figureFont...
                        ,'Interpreter','TeX'...
                        ,'Fontweight','bold'...
                         );
        hold off
        axis([xl yl]);
    end
    a = colorbar('Position', [0.92 0.11 0.02 .81]...
                 ,'Limits',[cMapMin,cMapMax]...
                 );

    cMapMid = mean([cMapMin,cMapMax]);
    if min(abs(cMapMid-Cs))<.002
        cMapMid = [];
    end
    newTicks = [cMapMin; cMapMax; cMapMid; Cs];
    [newTicks, Idx] = sort(newTicks);
    % a.Ticks = newTicks;
    webCodes = webCodes';
    if isempty(cMapMid)
        newTickLabels = [{num2str(cMapMin);
                          num2str(cMapMax)};
                          webCodes(:)];  
    else

        newTickLabels = [{num2str(cMapMin);
                          num2str(cMapMax);
                          num2str(cMapMid)};
                          webCodes(:)];
    end

    % a.TickLabels = newTickLabels(Idx);
    set(a, 'ytick', newTicks);
    set(a, 'yticklabel', newTickLabels(Idx));
    caxis([cMapMin,cMapMax])
    title(a,'C')

    arrayfun(@(x) set(x,'FontName',figureFont),get(localFigFirst,'children'));
    set(localFigFirst, 'paperunits', 'inches'...
                     , 'papertype','<custom>'...
                     , 'papersize', [img_width, img_height]...
                     , 'paperposition', [0, 0, img_width, img_height]...
        )
    fname = sprintf('../figures/initialProps%sLinkage_test.png',linkageType);
    print (fname) -r600 -tight
    
end

if plotsToPlot(2)||plotsToPlot(3)||plotsToPlot(4)
    allNodeLevels = propsGlobal(:,globalCol('S'));
    endings = [find(diff(allNodeLevels)>0); 
                numel(allNodeLevels)-1];
    beginnings = [1; endings(1:end-1)+1];
    endings = endings -5;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw plots of the agglomerated webs, size on x-axis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(2)    
      AR = 9/5.5;
    img_width = 18;
    img_height = img_width/AR;
    # Setup geometry of figure:
    left_offset = 0.5;
    bot_offset = 0.04;
    h_fig_sep = .02;
    v_fig_sep = 0.02;
    v_rem = 1 - bot_offset - 3*v_fig_sep;
    h_rem = 1 - .11 - h_fig_sep*3 - left_offset;
    fig_width = h_rem/4;
    fig_heighth = v_rem/3;
    localFigAllSizes = figure('units','inches'...
                           ,'position',[0,0,img_width,img_height]); ;
                           
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
        freesThisProperty = localMeans(:,localMeansCol(localPropsPlottedFree{jj}));
        parasThisProperty = localMeans(:,localMeansCol(localPropsPlottedPara{jj}));
        diffsThisProperty = freesThisProperty-parasThisProperty;
        [meanF_jj] = fake_grpstats(freesThisProperty,allNodeLevels,{@(x) mean(x)});
        [meanP_jj] = fake_grpstats(parasThisProperty,allNodeLevels,{@(x) mean(x)});
        meanFAll(jj,:) = meanF_jj;
        meanPAll(jj,:) = meanP_jj;
        [meanD_jj, numelD_jj, stdD_jj] = ...
        fake_grpstats(diffsThisProperty, allNodeLevels...
                    ,{   @(x) nanmean(x)...
                        ,@(x) sum(isfinite(x))...
                        ,@(x) nanstd(x)...
                     }...
                );
        #[meanD_jj,numelD_jj,stdD_jj] = grpstats(diffsThisProperty,allNodeLevels,{@mean,@numel,@std});
        ts = meanD_jj./(stdD_jj./sqrt(numelD_jj));
        ts(numelD_jj==1) = nan;
        Ps = (1-tcdf(abs(ts),numelD_jj-1))*2;

        allTs(:,jj) = ts;
        allPs(:,jj) = Ps;
        subplot(3,4,jj);

        hAll{jj} = plot(xAllLevels,meanP_jj - meanF_jj);
        set(hAll{jj},{'LineWidth'},{5});

        for ii = 1:nWebs
            hold on
            hEveryWeb{jj,ii} = plot(allNodeLevels(beginnings(ii):endings(ii)),parasThisProperty(beginnings(ii):endings(ii)) - freesThisProperty(beginnings(ii):endings(ii))...
                );
            set(hEveryWeb{jj,ii},{'LineWidth'},{.5})
            set(hEveryWeb{jj,ii},{'Color'},{[.7,.7,.7]})
            hold off

        end

        title(fancyLocalNames{jj},'Interpreter','tex')
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
        # This is the grey/pink shading.
        h = plot([xAllLevels';xAllLevels'],[zeros(size(meanPAll(jj,:)));meanPAll(jj,:) - meanFAll(jj,:)]);
        set(h, {'color'}, colsAllLevels(:,jj));
        set(h,'linewidth',3)
        #uistack(h,'bottom')
        #cellfun(@(x)uistack(x,'top'),hEveryWeb(jj,:),'UniformOutput',false);
        #uistack(hAll{jj},'top')
        #rl = refline(0,1);
        #rl.Color = [.5,.5,.5];
        #rl.LineStyle = '--';
        yl = ylim();
        plot([Ss(end);Ss(end)],[yl(1) ;yl(2)],'Color',[.6,.6,.6])
        axis([0,max(Ss) ylim]);
        ylim(yl);
        hold off
        
    end
    
    # arrayfun(@(x) set(x,'FontName',figureFont),localFigAllSizes.Children);
    set(localFigAllSizes, 'Units', 'Inches', 'Position', [0, 0, 17, 10])
    print ../figures/allPropsBySizeMaxLinkage.png -dpng -r0
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This sets up some important distance stuff.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(3)||plotsToPlot(4)||plotsToPlot(5)
    allDistances = localMeans(:,localMeansCol('minDistance'));
    nEdges = 81;
    nBins = nEdges-1;
    dMin = 0;
    dMax = 0.8;
    distancesToPlot = linspace(dMin,dMax,nEdges);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw plots of the agglomerated webs, minimum node distance on x-axis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figureFont =  'Times New Roman';
if plotsToPlot(3)
    AR = 9/5.5;
    img_width = 18;
    img_height = img_width/AR;
    % Setup geometry of figure:
    left_offset = 0.05;
    bot_offset = 0.04;
    h_fig_sep = .02;
    v_fig_sep = 0.02;
    v_rem = 1 - bot_offset - 3*v_fig_sep;
    h_rem = 1 - .11 - h_fig_sep*3 - left_offset;
    fig_width = h_rem/4;
    fig_heighth = v_rem/3;
    
    localFigDistance = figure('units','inches'...
                           ,'position',[0,0,img_width,img_height]); 
                      
    hDistances = cell(nLocalPropsPlotted,nWebs);
    hDistanceMean = cell(nLocalPropsPlotted,1);
    plotByMinDistance = nan(nEdges*nWebs,numel(localPropsPlotted));
    distanceGroups = repmat((1:nEdges)',nWebs,1);
    
    for ii = 1:nWebs
        distances_ii = allDistances(beginnings(ii):endings(ii));

        valuesFree_ii = localMeans(beginnings(ii):endings(ii),freeIndices);
        valuesPara_ii = localMeans(beginnings(ii):endings(ii),paraIndices); 
        diffs_ii = valuesPara_ii - valuesFree_ii;

        lastDistanceGreaterThan = sum(distances_ii<=distancesToPlot)';
        diffsToPlot_ii = diffs_ii(lastDistanceGreaterThan,:);

        plotByMinDistance(((1:nEdges)+nEdges*(ii-1)),:) = diffsToPlot_ii;

    end
    
    [groupedDistanceMeans, groupedDistanceStd, groupedDistanceN] = ...
        fake_grpstats(plotByMinDistance, distanceGroups...
                    ,{   @(x) nanmean(x)...
                        ,@(x) nanstd(x)...
                        ,@(x) sum(isfinite(x))...
                     }...
                );
            
    groupedDistanceTs = groupedDistanceMeans./(groupedDistanceStd./sqrt(groupedDistanceN));
    groupedDistancePs = 2-2*tcdf(abs(groupedDistanceTs),groupedDistanceN-1);
    
    #save(sprintf('../../Chapter3/code/groupedDistances%sLinkge.mat','max'),...
    #    'groupedDistanceMeans','groupedDistanceStd','groupedDistanceN'...
    #   ,'groupedDistanceTs','groupedDistancePs','plotByMinDistance'...
    #   ,'distancesToPlot');
   
   for jj = 1:nLocalPropsPlotted
   subplot(3,4,jj);
    hold on
        for ii = 1:nWebs
            thisWeb = (1:nEdges)+(ii-1)*nEdges;
            plotByMinDistanceThisWeb = plotByMinDistance(thisWeb,jj);
            hDistances{jj,ii} = plot(distancesToPlot,plotByMinDistanceThisWeb...
                ,'LineWidth',.5...
                ,'Color',[.5,.5,.5]);
        end
    hDistanceMean{jj} = plot(distancesToPlot,groupedDistanceMeans(:,jj)...
            ,'LineWidth', 2 ...
            ,'Color',[0 0 1]);
    hold off
    title(fancyLocalNames{jj},'Interpreter','tex')
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
        hDist = plot([distancesToPlot(2:end);distancesToPlot(2:end)],[zeros(size(distancesToPlot(2:end))); groupedDistanceMeans(2:end,jj)']);
        set(hDist, {'color'}, colsAllBins(2:end,jj));
        set(hDist,'linewidth',3
                 ,'clipping','on')
        #uistack(hDistanceMean,'top')
        #uistack(hDist,'bottom')
        #rl = refline(0,0);
        #rl.Color = [.5,.5,.5];
        #rl.LineStyle = '--';
        rl = line(distancesToPlot,zeros(size(distancesToPlot)));
        set(rl,'Color',[0.5, 0.5, 0.5]...
              ,'LineStyle','--')
        yl = ylim();
        ylim(yl);
        if jj > 8
            xlabel('Min. Cluster Distance')
        endif
        
        if mod(jj,4) == 1
            ylabel('Parasite - Carnivore','Interpreter','tex')
        endif
        hold off

    end
    cellfun(@(x) copyobj(x, get(x,'parent')), hDistances);
    cellfun(@(x) copyobj(x, get(x,'parent')), hDistanceMean); 
    arrayfun(@(x) set(x,'FontName',figureFont ...
                       , 'FontSize',14),get(localFigDistance,'children'));
    set(localFigDistance, 'paperunits', 'inches'...
                        , 'papertype','<custom>'...
                        , 'papersize', [img_width, img_height]...
                        , 'paperposition', [0, 0, img_width, img_height]...
        )
    fname = sprintf('../figures/distanceProps.png');
    print (fname) -r600 -tight
    
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
        title(globalPropsPlotted{jj},'Interpreter','tex')
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
        carnCarnJac = cellfun(@(s,t) triu(s(t>=0.5,t>=0.5),1),usedProps{ii}.jac(:),usedProps{ii}.carn(:),'UniformOutput',fpalse);
        
            
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Draw plots of the initial Species webs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(6)
    # Load the variables: OG_global_props & OG_local_means & OG_local_props
    load -v7 og_local_means.mat
    
    AR = 9/5.5;
    img_width = 18;
    img_height = img_width/AR;
    # Setup geometry of figure:
    left_offset = 0.05;
    bot_offset = 0.04;
    h_fig_sep = .02;
    v_fig_sep = 0.02;
    v_rem = 1 - bot_offset - 3*v_fig_sep;
    h_rem = 1 - .11 - h_fig_sep*3 - left_offset;
    fig_width = h_rem/4;
    fig_heighth = v_rem/3;
    
    localFigFirst = figure('units','inches'...
                           ,'position',[0,0,img_width,img_height]); 
                           
    colormap(summer);
    cm_ = get(localFigFirst,'colormap');
    set(localFigFirst,'colormap',cm_(size(cm_,1):-1:1,:));
    
    
    Ss = OG_global_props(:,globalCol('S')-1);
    Cs = OG_global_props(:,globalCol('C')-1);
    # selectThisLevelLocal = arrayfun(@(x,y) (propsLocal(:,localCol('S'))==y)&(propsLocal(:,localCol('web'))==x),webs',Ss,'UniformOutput',false);
    # initialWebs = sum(cell2mat(selectThisLevelLocal'),2)>0;
    # initialWebsPropsLocal = propsLocal(initialWebs,:);

    # initCarn = (initialWebsPropsLocal(:,localCol('free'))>0.5);
    # initPara = (initialWebsPropsLocal(:,localCol('para'))>=0.5);
    # initCarnOrPara = initCarn|initPara;

    cMapMin = min(OG_global_props(:,globalCol('C')-1))*.99-.001;
    cMapMax = max(OG_global_props(:,globalCol('C')-1))*1.01+.001;
    colMatrix = summer(100);    
    colMatrix(size(colMatrix,1):-1:1,:);
    cVals = linspace(cMapMin,cMapMax,100);

    
    tDiffs = zeros(nLocalPropsPlotted,1);
    pDiffs = zeros(nLocalPropsPlotted,1);
    xMeans = zeros(nLocalPropsPlotted,1);
    yMeans = zeros(nLocalPropsPlotted,1);
    
    rC = zeros(nLocalPropsPlotted,1);
    rS = zeros(nLocalPropsPlotted,1);
    
    for jj = 1:nLocalPropsPlotted
        figure(localFigFirst);
        xs = OG_local_means(:,localMeansCol(localPropsPlottedFree{jj})-1);
        ys = OG_local_means(:,localMeansCol(localPropsPlottedPara{jj})-1);
        gps = 1:6; #OG_local_means(:,localMeansCol('web'));
        markSizes = sqrt(Ss-100)*30;
        markColors = colMatrix(arrayfun(@(x)find(x<cVals,1),Cs),:);
        fig_col = mod(jj-1,4) + 1;
        fig_row = floor((jj-1)/4)+1;
        h_pos = left_offset + (fig_col-1)*(fig_width + h_fig_sep)
        v_pos = bot_offset + (3 - fig_row)*(fig_heighth + v_fig_sep)
        subplot(3,4,jj,'position',[h_pos, v_pos, fig_width, fig_heighth]);
        h = scatter(xs,ys,markSizes,markColors,'filled');
        rC(jj) = corr(ys-xs,Cs);
        rS(jj) = corr(ys-xs,Ss);
        if jj > 8
            xlabel('Carnivore Average')
        else
            xlabel('')
        end
        if mod(jj,4) == 1
            ylabel('Parasite Average')
        else
            ylabel('')
        end
        title(fancyLocalNames{jj}...
            ,'Interpreter','tex'...
            ,'FontSize',14)
        xl = xlim;
        yl = ylim;
        xMeans(jj) = mean(xs);
        yMeans(jj) = mean(ys);
        diffs = ys-xs;
        meanDiff = mean(diffs);
        stdDiff = std(diffs)/sqrt(numel(diffs));
        tDiff = meanDiff./stdDiff;
        tDiffs(jj) = tDiff;
        pDiff = 1-tcdf(abs(tDiff),numel(diffs)-1);
        pDiffs(jj) = pDiff;
        lb = min(xl(1),yl(1));
        ub = max(xl(2),yl(2));
        
        %rl = refline(1,0);
        %rl.Color = textColFailed;
        % uistack(h,'top')
        hold on

        rl = line([lb ub],[lb ub],'color',textColFailed);
        h = scatter(xs,ys,markSizes,markColors);
        axis([lb ub lb ub]);
    end
    
    [pDiffsSorted,Idx] = sort(pDiffs);
    [~,IdxIdx] = sort(Idx);
    rejected = false(nLocalPropsPlotted,1);
    rejected(Idx) = pDiffsSorted<alphaSeq';
    
    tC = rC*2./sqrt(1-rC.^2);
    tS = rS*2./sqrt(1-rS.^2);
    pC = (1-tcdf(abs(tC),4))*2;
    pS = (1-tcdf(abs(tS),4))*2;
    
    CS_tab = [rC tC pC rS tS pS];
    data_cell = horzcat(fancyLocalNames',mat2cell(CS_tab,repmat([1],12,1),repmat([1],6,1)));
    fid = fopen('interaction_cs.csv','w')
    fprintf(fid,',corr(C),t(C),P(C),corr(S),t(S),P(S)\n')
    fprintf(fid,'%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n',data_cell'{:})
    fclose(fid)
    for jj = 1:nLocalPropsPlotted
       fig_col = mod(jj-1,4) + 1;
       fig_row = floor((jj-1)/4)+1;
       h_pos = left_offset + (fig_col-1)*(fig_width + h_fig_sep);
       v_pos = bot_offset + (3 - fig_row)*(fig_heighth + v_fig_sep);
       subplot(3,4,jj)#,'position',[h_pos, v_pos, fig_width, fig_heighth]);
        hold on
        xl = xlim();
        yl = ylim();
        
        if rejected(jj)
            textCol = textColRejected;
            plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
            plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
            meanLabel = sprintf('t=%.3g\nP=%.3g<%.3g',tDiffs(jj),pDiffs(jj),alphaSeq(IdxIdx(jj)));
        else
            textCol = textColFailed;
            plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
            plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
            meanLabel = sprintf('t=%.3g\nP=%.3g\\geq%.3g',tDiffs(jj),pDiffs(jj),alphaSeq(IdxIdx(jj)));
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
                        ,'FontSize',14 ...
                        ,thisAlign{:}...
                        ,'Margin',1 ...
                        ,'backgroundColor','white'...
                        ,'FontName',figureFont...
                        ,'Interpreter','TeX'...
                        ,'Fontweight','bold'...
                         );
        hold off
        axis([xl yl]);
    end
    a = colorbar('Position', [0.92 0.11 0.02 .81]...
                 ,'Limits',[cMapMin,cMapMax]...
                 );

    cMapMid = mean([cMapMin,cMapMax]);
    if min(abs(cMapMid-Cs))<.002
        cMapMid = [];
    end
    newTicks = [cMapMin; cMapMax; cMapMid; Cs];
    [newTicks, Idx] = sort(newTicks);
    % a.Ticks = newTicks;
    webCodes = webCodes';
    if isempty(cMapMid)
        newTickLabels = [{num2str(cMapMin);
                          num2str(cMapMax)};
                          webCodes(:)];  
    else

        newTickLabels = [{num2str(cMapMin);
                          num2str(cMapMax);
                          num2str(cMapMid)};
                          webCodes(:)];
    end

    % a.TickLabels = newTickLabels(Idx);
    set(a, 'ytick', newTicks);
    set(a, 'yticklabel', newTickLabels(Idx));
    caxis([cMapMin,cMapMax])
    title(a,'C')

    arrayfun(@(x) set(x,'FontName',figureFont),get(localFigFirst,'children'));
    set(localFigFirst, 'paperunits', 'inches'...
                     , 'papertype','<custom>'...
                     , 'papersize', [img_width, img_height]...
                     , 'paperposition', [0, 0, img_width, img_height]...
        )
    fname = sprintf('../figures/initial_species_Props%sLinkage_test.png',linkageType);
    print (fname) -r600 -tight
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Histograms of the initial Species webs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(8)
    # Load the variables: global_props & local_mean_props & local_props
    load -v7 og_local_means.mat
    
    AR = 9/5.5;
    img_width = 18;
    img_height = img_width/AR;
    # Setup geometry of figure:
    left_offset = 0.05;
    bot_offset = 0.04;
    h_fig_sep = .02;
    v_fig_sep = 0.02;
    v_rem = 1 - bot_offset - 3*v_fig_sep;
    h_rem = 1 - 2*left_offset - h_fig_sep*3;
    fig_width = h_rem/4;
    fig_heighth = v_rem/3;
    
    localFigFirst = figure('units','inches'...
                           ,'position',[0,0,img_width,img_height]); 
                           
    colormap(summer);
    cm_ = get(localFigFirst,'colormap');
    set(localFigFirst,'colormap',cm_(size(cm_,1):-1:1,:));
    
    
    Ss = OG_global_props(:,globalCol('S')-1);
    Cs = OG_global_props(:,globalCol('C')-1);
    # selectThisLevelLocal = arrayfun(@(x,y) (propsLocal(:,localCol('S'))==y)&(propsLocal(:,localCol('web'))==x),webs',Ss,'UniformOutput',false);
    # initialWebs = sum(cell2mat(selectThisLevelLocal'),2)>0;
    # initialWebsPropsLocal = propsLocal(initialWebs,:);

    # initCarn = (initialWebsPropsLocal(:,localCol('free'))>0.5);
    # initPara = (initialWebsPropsLocal(:,localCol('para'))>=0.5);
    # initCarnOrPara = initCarn|initPara;

    # cMapMin = min(global_props(:,globalCol('C')-1))*.99-.001;
    # cMapMax = max(global_props(:,globalCol('C')-1))*1.01+.001;
    # colMatrix = summer(100);    
    # colMatrix(size(colMatrix,1):-1:1,:);
    # cVals = linspace(cMapMin,cMapMax,100);

    
    tDiffs = zeros(nLocalPropsPlotted,1);
    pDiffs = zeros(nLocalPropsPlotted,1);
    xMeans = zeros(nLocalPropsPlotted,1);
    yMeans = zeros(nLocalPropsPlotted,1);
    
    rC = zeros(nLocalPropsPlotted,1);
    rS = zeros(nLocalPropsPlotted,1);
    OG_local_props(isnan(OG_local_props)) = -1;
    para = OG_local_props(:,localCol('para') - 1);
    carn = OG_local_props(:,localCol('free') - 1);
    var_hists = cell(nLocalPropsPlotted,2);
    for jj = 1:nLocalPropsPlotted
        figure(localFigFirst);
        X = OG_local_props(:,localCol(localPropsPlotted{jj})-1);
        X_par = X(para>0);
        X_carn = X(carn>0);
        fig_col = mod(jj-1,4) + 1;
        fig_row = floor((jj-1)/4)+1;
        h_pos = left_offset + (fig_col-1)*(fig_width + h_fig_sep)
        v_pos = bot_offset + (3 - fig_row)*(fig_heighth + v_fig_sep)
        sp = subplot(3,4,jj,'position',[h_pos, v_pos, fig_width, fig_heighth]);
        hold on
        [~, bins] = hist(X);
        bin_width = diff(bins,15);
        [n1, x1] = hist(X_par, bins, 'facecolor','r'...
                             , 'edgecolor','k');
        [n2, x2] = hist(X_carn,bins, 'facecolor','b'...
                   ,'edgecolor','k');
        xx = [x1; x1];
        xx = xx(2:end) - bin_width(1)/2;
        nn1 = [n1; n1];
        nn1 = nn1(1:end-1);
        nn2 = [n2; n2];
        nn2 = nn2(1:end-1);
        xx = [xx(1) xx xx(end)];
        nn1 = [0 nn1 0];
        nn2 = [0 nn2 0];
        h1=fill(xx,nn1,"red");
        h2=fill(xx,nn2,"blue");
        xlim([min(xx), max(xx)]);
        set(h1,'facealpha',0.5)
        set(h2,'facealpha',0.5)
        
        if jj > 8
            xlabel('Property Value')
        else
            xlabel('')
        end
        if mod(jj,4) == 1
            ylabel('Frequency')
        else
            ylabel('')
        end
        title(fancyLocalNames{jj}...
            ,'Interpreter','tex'...
            ,'FontSize',14);
    end
    legend('Parasites', 'Carnivores')

    arrayfun(@(x) set(x,'FontName',figureFont),get(localFigFirst,'children'));
    set(localFigFirst, 'paperunits', 'inches'...
                     , 'papertype','<custom>'...
                     , 'papersize', [img_width, img_height]...
                     , 'paperposition', [0, 0, img_width, img_height]...
        )
    fname = sprintf('../figures/initial_species_Hists_test.png',linkageType);
    print (fname) -r600 -tight
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Histograms of the initial Trophic webs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotsToPlot(8)
   
    AR = 9/5.5;
    img_width = 18;
    img_height = img_width/AR;
    # Setup geometry of figure:
    left_offset = 0.05;
    bot_offset = 0.04;
    h_fig_sep = .02;
    v_fig_sep = 0.02;
    v_rem = 1 - bot_offset - 3*v_fig_sep;
    h_rem = 1 - 2*left_offset - h_fig_sep*3;
    fig_width = h_rem/4;
    fig_heighth = v_rem/3;
    
    localFigFirst = figure('units','inches'...
                           ,'position',[0,0,img_width,img_height]); 
                           
    Ss = OG_global_props(:,globalCol('S')-1);
    # selectThisLevelLocal = arrayfun(@(x,y) (propsLocal(:,localCol('S'))==y)&(propsLocal(:,localCol('web'))==x),webs',Ss,'UniformOutput',false);
    # initialWebs = sum(cell2mat(selectThisLevelLocal'),2)>0;
    # initialWebsPropsLocal = propsLocal(initialWebs,:);

    # initCarn = (initialWebsPropsLocal(:,localCol('free'))>0.5);
    # initPara = (initialWebsPropsLocal(:,localCol('para'))>=0.5);
    # initCarnOrPara = initCarn|initPara;

    # cMapMin = min(global_props(:,globalCol('C')-1))*.99-.001;
    # cMapMax = max(global_props(:,globalCol('C')-1))*1.01+.001;
    # colMatrix = summer(100);    
    # colMatrix(size(colMatrix,1):-1:1,:);
    # cVals = linspace(cMapMin,cMapMax,100);

    initial_webs = propsLocal(:,localCol('counter'))==1;
    web1 = propsLocal(:,localCol('web')) == 1;
    para = propsLocal(initial_webs, localCol('para'));
    carn = propsLocal(initial_webs, localCol('free'));
    
    var_hists = cell(nLocalPropsPlotted,2);
    for jj = 1:nLocalPropsPlotted
        figure(localFigFirst);
        localPropsPlotted{jj}
        X = propsLocal(initial_webs,localCol(localPropsPlotted{jj}));
        X(isnan(X)) = -1;
        X_par = X(para>0);
        X_carn = X(carn>0);
        mean(X_par)
        mean(X_carn)
        fig_col = mod(jj-1,4) + 1;
        fig_row = floor((jj-1)/4)+1;
        h_pos = left_offset + (fig_col-1)*(fig_width + h_fig_sep);
        v_pos = bot_offset + (3 - fig_row)*(fig_heighth + v_fig_sep);
        sp = subplot(3,4,jj,'position',[h_pos, v_pos, fig_width, fig_heighth]);
        hold on
        [~, bins] = hist(X,15);
        bin_width = diff(bins);
        [n1, x1] = hist(X_par, bins, 'facecolor','r'...
                             , 'edgecolor','k');
        [n2, x2] = hist(X_carn,bins, 'facecolor','b'...
                   ,'edgecolor','k');
        xx = [x1; x1];
        xx = xx(2:end) - bin_width(1)/2;
        nn1 = [n1; n1];
        nn1 = nn1(1:end-1);
        nn2 = [n2; n2];
        nn2 = nn2(1:end-1);
        xx = [xx(1) xx xx(end)];
        nn1 = [0 nn1 0];
        nn2 = [0 nn2 0];
        h1=fill(xx,nn1,"red");
        h2=fill(xx,nn2,"blue");
        xlim([min(xx), max(xx)]);
        
        set(h1,'facealpha',0.5)
        set(h2,'facealpha',0.5)
        
        if jj > 8
            xlabel('Property Value')
        else
            xlabel('')
        end
        if mod(jj,4) == 1
            ylabel('Frequency')
        else
            ylabel('')
        end
        title(fancyLocalNames{jj}...
            ,'Interpreter','tex'...
            ,'FontSize',14);
    end
    legend('Parasites', 'Carnivores')

    arrayfun(@(x) set(x,'FontName',figureFont),get(localFigFirst,'children'));
    set(localFigFirst, 'paperunits', 'inches'...
                     , 'papertype','<custom>'...
                     , 'papersize', [img_width, img_height]...
                     , 'paperposition', [0, 0, img_width, img_height]...
        )
    fname = sprintf('../figures/initial_trophic_Hists_test.png',linkageType);
    print (fname) -r600 -tight
    
end
localMeans(localMeans(:,localMeansCol('counter'))==1,:)
web1_props = propsLocal(web1 & initial_webs,:);