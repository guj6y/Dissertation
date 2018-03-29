%Biggest question here is HOW DO I REACH THESE KEEDS???

%jk, biggest question is how to align the webs so that we can get a big fat
%set of species to look at...
%option 1: By web size
%   Pros: easy-peasy, somewhat standard control
%   Cons: agglomerated webs must be different from normal webs, right?;
%         
%option 2: By minimum node distance
%   Pros: Agrees with the process; accounts for differences in the very
%   process i am testing...
%   Cons: Harder to line things up right
%         
%Another issue is whether to compare 
% 1. parasites to free-livers, 
% 2. parasites to carnivores,
% 3. either of the above in an invertebrate-basal web
% For whatever reason, I'm not really a big fan of the invertebrate-basal
% web so I'm just comparing to carnivores since option 1 is pretty clearly
% wrong anyway.

%SPOILER ALERT: initial webs should be pretty opaque by these measures; the
%**hope** is that later in the agglomeration sequences a pattern deus ex
%machina's itself and slaps me across the face

%Linear Regression
%Logit model
%Classification Trees
%knn?

%Setting up.
clear

for linkageTypes = {'Max','Min'}
linkageType = linkageTypes{:};
load('data/Processed/webGeneration.mat')
T = load(sprintf('AgglomerationProps%sLinkage.mat',linkageType));

propsLocal = T.propsLocal;
localCol = T.localCol;
propsGlobal = T.propsGlobal;
globalCol = T.globalCol;
localMeans = T.propsLocalMeans;
localMeansCol = T.localMeansCol;
clear T
%In case you were wondering.
%{
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
                ,'free'       ...12 actually, carn
                ,'S'          ...13
                ,'C'          ...14
                ,'size'       ...15
                ,'ccCyc'      ...16
                ,'ccMid'      ...17
                ,'ccIn'       ...18
                ,'ccOut'      ...19
                };
%}
variableNames = {'vul'        ...1      2
                ,'gen'        ...2      3
                ,'patl'       ...3      4
                ,'pr'         ...4      6
                ,'btwn'       ...5      7
                ,'ecoBtwn'    ...6      8
                ,'meanVulPrey'...7      9
                ,'meanGenPred'...8      10
                ,'ccCyc'      ...9      15
                ,'ccMid'      ...10     16
                ,'ccIn'       ...11     17
                ,'ccOut'      ...12     18
                ...,'meanGenPred2'...8     19
                ...,'ccCyc2'      ...9     20
                ...,'ccMid2'      ...10    21
                ...,'ccIn2'       ...11    22
                ...,'ccOut2'      ...12    23
                ,'para'       ...13
                };
set(groot, 'defaultAxesTickLabelInterpreter','LaTeX'); 
set(groot, 'defaultLegendInterpreter','LaTeX');

fancyNames = {'$v$'...                  1
             ,'$g$'...                  2
             ,'$TL$'...                 3
             ,'$\lambda$'...            4
             ,'$C_B$'...                5
             ,'$C_{EB}$'...             6
             ,'$v_r$'...                7
             ,'$g_c$'...                8
             ,'$\gamma^{cyc}$'...       9
             ,'$\gamma^{mid}$'...       10
             ,'$\gamma^{snk}$'...       11
             ,'$\gamma^{src}$'...       12
             };
nPred = numel(variableNames);


%propsThatCanBeNan = propsLocal(:,[10, 15,16,17,18]);
%obsThatAreNan = isnan(propsThatCanBeNan);
%propsThatCanBeNan(obsThatAreNan) = 0;
%propsLocal(:,[10, 15,16,17,18]) = propsThatCanBeNan;
%propsLocal = [propsLocal obsThatAreNan];

figureFont = 'CMU Serif';

sMin = min(propsLocal(:,localCol('S')));
sMax = max(propsLocal(:,localCol('S')));

carnAndPara = (propsLocal(:,localCol('free'))>0.5)|(propsLocal(:,localCol('para'))>=0.5);
webs = 1:6;
selectThisLevelGlobal = arrayfun(@(x) find(propsGlobal(:,globalCol('web'))==x,1),webs);
Ss = propsGlobal(selectThisLevelGlobal,globalCol('S'));
selectThisLevelLocal = arrayfun(@(x,y) (propsLocal(:,localCol('S'))==y)&(propsLocal(:,localCol('web'))==x),webs',Ss,'UniformOutput',false);
%
%
initialWebs = sum(cell2mat(selectThisLevelLocal'),2)>0;
initialWebsPropsLocal = propsLocal(initialWebs,:);
initCarn = (initialWebsPropsLocal(:,localCol('free'))>0.5);
initPara = (initialWebsPropsLocal(:,localCol('para'))>=0.5);
initCarnOrPara = initCarn|initPara;
%
%
X0 = propsLocal(:,[2:4 6:10, 16:19]);
X0(isnan(X0)) = 0;
Y0 = propsLocal(:,localCol('para'))>=0.5; %categorical
%
XForBig = initialWebsPropsLocal(initCarnOrPara,[2:4 6:10, 16:19]);
XForBig(isnan(XForBig)) = 0;
YForBig = initialWebsPropsLocal(initCarnOrPara,localCol('para'))>=0.5;
fParaBig = mean(YForBig);

nBig = numel(YForBig);
ogNodeList = 1:nBig;
nBigTest = 100;
classificationPlot = figure;

impVars = [1:12];
%impVars = [3 5 8 10 12];
XImp = XForBig(:,impVars);
varNames = variableNames(impVars);
fancyVarNames = fancyNames(impVars);
impVars = [impVars 0];
nImpVars = numel(impVars);

relErrorMeans = zeros(nImpVars-1,2);
relErrorStds = zeros(nImpVars-1,2);
truePositiveMeans = zeros(nImpVars-1,2);
truePositiveStds = zeros(nImpVars-1,2);
bigFatTree = fitctree(XImp,YForBig...
                        ,'PredictorNames',varNames...
                        ...,'Prune','on'...
                        ,'MaxNumSplits',4 ...
                        ...,MinleafSize',10...
                    );
bigFatFracPara = mean(YForBig);

for  ii= 1:nImpVars

    knockMeDown = impVars(ii);
    if knockMeDown >0
    distbnOfPredictors = zeros(nPred-1,1);
    %Test how it does without the big 5 predictors: 3,5,8,10,12:
    %
    XOut = XImp;
    XIn = XOut(:,ii);
    XOut(:,ii) = [];

    bigTreeErrors = zeros(nBigTest,2);
    bigTreeTruePositive = zeros(nBigTest,2);
    bigTreeRelErrors = zeros(nBigTest,2);
    varNamesOut = varNames;
    varNamesIn = varNamesOut(ii);
    varNamesOut(ii) = [];
    else
        fullError = zeros(nBigTest,1);
        fullTruePositive = zeros(nBigTest,1);
        fullRelError = zeros(nBigTest,1);
    end
    
    %}
    for kk = 1:nBigTest
            kk

            trainingIndices = randsample(ogNodeList,round(0.5*nBig));
            pickTrain = false(nBig,1);
            pickTrain(trainingIndices) = true;
            pickTest = ~pickTrain;
            if knockMeDown >0
                xTrainOut = XOut(pickTrain,:);
                xTrainIn = XIn(pickTrain);
                yTrain = YForBig(pickTrain,:);

                xTestOut = XOut(pickTest,:);
                xTestIn = XIn(pickTest);
                yTest = YForBig(pickTest,:);

                bigTreeIn = fitctree(xTrainIn,yTrain...
                        ,'PredictorNames',varNamesIn...
                        ...,'Prune','on'...
                        ,'MaxNumSplits',4 ...
                        ...,MinleafSize',10...
                    );

                 bigTreeOut = fitctree(xTrainOut,yTrain...
                        ,'PredictorNames',varNamesOut...
                        ...,'Prune','on'...
                        ,'MaxNumSplits',4 ...
                        ...,MinleafSize',10...
                    );

                fracParaTrain = mean(yTrain);

                usedInTree = bigTreeIn.CutPredictor;
                
                knownFractionErrors = 2*fracParaTrain*(1-fracParaTrain);

                bigTreeErrors(kk,1) = (1 - sum(bigTreeIn.predict(xTestIn)==yTest)/numel(yTest));
                bigTreeRelErrors(kk,1) = (1 - sum(bigTreeIn.predict(xTestIn)==yTest)/numel(yTest))/knownFractionErrors;
                bigTreeTruePositive(kk,1) = 1-sum(bigTreeIn.predict(xTestIn)&yTest)/sum(yTest);

                bigTreeErrors(kk,2) = (1 - sum(bigTreeOut.predict(xTestOut)==yTest)/numel(yTest));
                bigTreeRelErrors(kk,2) = (1 - sum(bigTreeOut.predict(xTestOut)==yTest)/numel(yTest))/knownFractionErrors;
                bigTreeTruePositive(kk,2) = 1-sum(bigTreeOut.predict(xTestOut)&yTest)/sum(yTest);
            else
                xTrain = XImp(pickTrain,:);
                yTrain = YForBig(pickTrain,:);

                xTest = XImp(pickTest,:);
                yTest = YForBig(pickTest,:);

                 bigTree = fitctree(xTrain,yTrain...
                        ,'PredictorNames',varNames'...
                        ...,'Prune','on'...
                        ,'MaxNumSplits',4 ...
                        ...,MinleafSize',10...
                    );

                fracParaTrain = mean(yTrain);

                
                knownFractionErrors = 2*fracParaTrain*(1-fracParaTrain);

                fullError(kk) = (1 - sum(bigTree.predict(xTest)==yTest)/numel(yTest));
                fullRelError(kk) = (1 - sum(bigTree.predict(xTest)==yTest)/numel(yTest))/knownFractionErrors;
                fullTruePositive(kk) = 1-sum(bigTree.predict(xTest)&yTest)/sum(yTest);
                
            end
    end
    if knockMeDown == 0
        fullErrorMean = mean(fullError);
        fullErrorStd = std(fullError);
        fullRelErrorMean = mean(fullRelError);
        fullRelErrorStd = std(fullRelError);
        fullTruePositiveMean = mean(fullTruePositive);
        fullTruePositiveStd = std(fullTruePositive);
    else
            
        relErrorMeans(ii,:) = mean(bigTreeRelErrors);
        relErrorStds(ii,:) = std(bigTreeRelErrors);
        truePositiveMeans(ii,:) = mean(bigTreeTruePositive);
        truePositiveStds(ii,:) = std(bigTreeTruePositive);

    
    end
end


classFigure = figure;
a1 = subplot(2,2,1);
plotABar(relErrorMeans(:,2),relErrorStds(:,2),fullRelErrorMean);
title('(a) Relative Error, Removed Variable')
xlabel('Variable Removed')
ylabel('Error Relative to Random Model')


a2 = subplot(2,2,2);
plotABar(relErrorMeans(:,1),relErrorStds(:,1),fullRelErrorMean);
title('(b) Relative Error, Isolated Variable')
xlabel('Variable Isolated')
ylabel('Error Relative to Random Model')


a3 = subplot(2,2,3);
plotABar(truePositiveMeans(:,2),truePositiveStds(:,2),fullTruePositiveMean);
title('(c) True Positive, Removed Variable')
xlabel('Variable Removed')
ylabel('Fraction of Parasites Misclassifid')

a4 = subplot(2,2,4);
plotABar(truePositiveMeans(:,1),truePositiveStds(:,1),fullTruePositiveMean);
title('(d) True Positive, Isolated Variable')
xlabel('Variable Isolated')
ylabel('Fraction of Parasites Misclassified')

arrayfun(@(x) set(x,'FontName','CMU Serif'),classFigure.Children);
%arrayfun(@(x) set(x,'Interpreter','LaTeX'),[a1.Children;a2.Children;a3.Children;a4.Children]);
arrayfun(@(x) set(x,'XTickLabels',fancyVarNames),classFigure.Children);
arrayfun(@(x) set(x,'XTickLabelRotation',20),classFigure.Children);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 13, 6.5]);
arrayfun(@(x) set(x,'XTick',1:12,'XLim',[0 13]),classFigure.Children)
impVars(end) = [];
print(sprintf('../figures/treeAnalysis%sLinkage.png',linkageType),'-dpng','-r0')
%}
%This is insanity. Hard enough to make one model; why try to make an
%arbitrary number of them. They were all worthless; the big one is 
%worseless on the agglomerated webs, too; so why did I even bother with
%that?
%They're not worthless now.. what happened again? The variables weren't
%properly normalized.


%Need to make an array of distances for EACH species; didn't include that
%and at this point adding it in will be a major pain in the ass.
allNodeLevels = propsGlobal(:,globalCol('S'));
allDistances = localMeans(:,localMeansCol('minDistance'));
allWebs = localMeans(:,localMeansCol('web'));
xAllLevels = unique(allNodeLevels);
endings = [find(diff(allNodeLevels)>0); numel(allNodeLevels)-1];
beginnings = [1; endings(1:end-1)+1];
endings = endings -5;
nEdges = 81;
nBins = nEdges-1;
dMin = 0;
dMax = 0.4;
distancesToPlot = linspace(dMin,dMax,nEdges);

plotByMinDistance = nan(nEdges*nWebs,12);
websByMinDistance = reshape(repmat(1:nWebs,nEdges,1),[],1);
distanceGroups = repmat((1:nEdges)',nWebs,1);
allMinDistances = zeros(nEdges,nWebs);

selectors = nan(nEdges*nWebs,4);


allSs = zeros(nEdges,nWebs);

for ii = 1:nWebs
    distances_ii = allDistances(beginnings(ii):endings(ii));
    nodeLevels_ii = allNodeLevels(beginnings(ii):endings(ii));
    findMinDistances = sum(distances_ii<=distancesToPlot);
    minDistances_ii = distances_ii(findMinDistances);
    minNodeLevels_ii = nodeLevels_ii(findMinDistances);
    
    %Find S that corresponds to each minimu distance.
    
    selectors((1:nEdges)+(ii-1)*nEdges,:) = ...
       [(1:nEdges)'  minNodeLevels_ii repmat(ii,nEdges,1) minDistances_ii];
    
end
save(sprintf('selectors%sLinkage.mat',linkageType),'selectors');
meanWebSizes = grpstats(selectors(:,2),selectors(:,1));

nAllSpecies = length(propsLocal);
nAllWebs = length(propsGlobal);
allDistances = nan(nAllSpecies,1);
theseDistanceGroups = zeros(nAllSpecies,1);



%This is actually pretty quick.
for ii = 1:nAllWebs
    
    %fprintf('%.1f%%\n',ii/nAllWebs*100);
    thisWeb = propsGlobal(ii,globalCol('web'));
    thisWebsS = propsGlobal(ii,globalCol('S'));
    thisWebsD = propsGlobal(ii,globalCol('d_J'));
    
    theseSpeciesForThisWeb = (propsLocal(:,localCol('web'))==thisWeb) & ...
        (propsLocal(:,localCol('S'))==thisWebsS);
    
    allDistances(theseSpeciesForThisWeb) = thisWebsD;
end

trees = cell(nEdges,1);


naddas  = nan(nEdges,1);


treeTruePositiveMean = nan(nEdges,1);
treeRelErrorsMean = nan(nEdges,1);
treeErrorsMean = nan(nEdges,1);

treeTruePositiveStd = nan(nEdges,1);
treeRelErrorsStd = nan(nEdges,1);
treeErrorsStd = nan(nEdges,1);

bigFatTreeTruePositiveMean = nan(nEdges,1);
bigFatTreeRelErrorsMean = nan(nEdges,1);
bigFatTreeErrorsMean = nan(nEdges,1);

bigFatTreeTruePositiveStd = nan(nEdges,1);
bigFatTreeRelErrorsStd = nan(nEdges,1);
bigFatTreeErrorsStd = nan(nEdges,1);

SbySpecies = propsLocal(:,localCol('S'));
webBySpecies = propsLocal(:,localCol('web'));

%predFrequency = zeros(nEdges,nPred-1);
shortLocalCol = containers.Map(variableNames,1:13);
nReps = round(logspace(2,3,81));
%nReps = round(linspace(sqrt(10),sqrt(100),81).^2);
try
    load(sprintf('predFrequency%sLinkage.mat',linkageType));
catch
    getPredTotal = 2000;
    predFrequency = zeros(getPredTotal,nPred-1);
    for ii = 1:nEdges
        sAndWeb = selectors(selectors(:,1) == ii,[2 3]);
        usesTheseSpecies = false(nAllSpecies,1);
        useTheseSpecies = carnAndPara&...
            sum((SbySpecies == sAndWeb(:,1)')&(webBySpecies==sAndWeb(:,2)'),2);
        varNames = variableNames(1:nPred-1);
        X = X0(useTheseSpecies,:);
        Y = Y0(useTheseSpecies);
        
        nTotal = length(Y);
        trainFraction = 0.8;
        nTrain = round(trainFraction*nTotal);
        nTest = nTotal-nTrain;
        
        for jj = 1:getPredTotal
            
            testThese = false(nTotal,1);
            testIndices = randsample(nTotal,nTest);
            testThese(testIndices) = true;
            xTest = X(testThese,:);
            yTest = Y(testThese,:);
            xTrain = X(~testThese,:);
            yTrain = Y(~testThese,:);
            
            tree = fitctree(xTrain,yTrain...
                    ,'PredictorNames',varNames...
                    ,'MaxNumSplits',4 ...
                        );
            fracParaTrain = mean(yTrain);
            
            predsUsed = tree.CutPredictor;
            predsUsed(cellfun(@isempty,predsUsed)) = [];
            predUsedIdx = cellfun(@(x) shortLocalCol(x),predsUsed);
            predFrequency(ii,predUsedIdx) = predFrequency(ii,predUsedIdx) + 1;
        end
    end
    save(sprintf('predFrequency%sLinkage.mat',linkageType));
    
end
    
predsUsedBin = false(nEdges,nPred-1);
trainFraction = 0.8;
for ii = 1:nEdges
    ii
    %Need to pick out the webs I'm using; useTheseWebs gives S, webNo for
    %all six webs.
    predFreqThisSize = predFrequency(ii,:);
    [predsSorted,useThesePreds] = sort(predFreqThisSize,'descend');
    useThesePreds = useThesePreds(1:4);
    %useThesePreds = 1:12;
    predsUsedBin(ii,useThesePreds) = true;
    sAndWeb = selectors(selectors(:,1) == ii,[2 3]);
    usesTheseSpecies = false(nAllSpecies,1);
    useTheseSpecies = carnAndPara&...
        sum((SbySpecies == sAndWeb(:,1)')&(webBySpecies==sAndWeb(:,2)'),2);
    varNames = variableNames(useThesePreds);
    X = X0(useTheseSpecies,useThesePreds);
    Y = Y0(useTheseSpecies);
    
    bigFatX0 = X0(useTheseSpecies,impVars);
    
    nTotal = length(Y);
    nTrain = round(trainFraction*nTotal);
    nTest = nTotal-nTrain;
    treeErrorsTemp = zeros(nReps(ii),1);
    treeRelErrorsTemp = zeros(nReps(ii),1);
    treeTruePositiveTemp = zeros(nReps(ii),1);
    
    bigFatTreeErrorsTemp = zeros(nReps(ii),1);
    bigFatTreeRelErrorsTemp = zeros(nReps(ii),1);
    bigFatTreeTruePositiveTemp = zeros(nReps(ii),1);
    
    for kk = 1:nReps(ii)
        testThese = false(nTotal,1);
        testIndices = randsample(nTotal,nTest);
        testThese(testIndices) = true;
        xTest = X(testThese,:);
        yTest = Y(testThese,:);
        xTrain = X(~testThese,:);
        yTrain = Y(~testThese,:);
        bigFatX = bigFatX0(testThese,:);
        
        tree = fitctree(xTrain,yTrain...
                    ,'PredictorNames',varNames...
                    ,'MaxNumSplits',4 ...
                        );
                    fracParaTrain = mean(yTrain);
        knownFractionErrors = 2*fracParaTrain*(1-fracParaTrain);
        treeErrorsTemp(kk) = (1 - sum(tree.predict(xTest)==yTest)/numel(yTest));
        treeRelErrorsTemp(kk) = (1 - sum(tree.predict(xTest)==yTest)/numel(yTest))/knownFractionErrors;
        treeTruePositiveTemp(kk) = 1-sum(tree.predict(xTest)&yTest)/sum(yTest);           
        
        bigFatTreeErrorsTemp(kk) = (1 - sum(bigFatTree.predict(bigFatX)==yTest/numel(yTest)));
        bigFatTreeRelErrorsTemp(kk) = (1 - sum(bigFatTree.predict(bigFatX)==yTest)/numel(yTest))/bigFatFracPara;
        bigFatTreeTruePositiveTemp(kk) = 1-sum(bigFatTree.predict(bigFatX)&yTest)/sum(yTest);           
    %{	
    predsUsed = tree.CutPredictor;
    predsUsed(cellfun(@isempty,predsUsed)) = [];
    predUsedIdx = cellfun(@(x) shortLocalCol(x),predsUsed);
    predFrequency(ii,predUsedIdx) = predFrequency(ii,predUsedIdx) + 1;
    %}
    end
    treeTruePositiveMean(ii) = mean(treeTruePositiveTemp);
    treeRelErrorsMean(ii) = mean(treeRelErrorsTemp);
    treeErrorsMean(ii) = mean(treeErrorsTemp);

    treeTruePositiveStd(ii) = std(treeTruePositiveTemp);
    treeRelErrorsStd(ii) = std(treeRelErrorsTemp);
    treeErrorsStd(ii) = std(treeErrorsTemp);
    
    bigFatTreeTruePositiveMean(ii) = mean(bigFatTreeTruePositiveTemp);
    bigFatTreeRelErrorsMean(ii) = mean(bigFatTreeRelErrorsTemp);
    bigFatTreeErrorsMean(ii) = mean(bigFatTreeErrorsTemp);

    bigFatTreeTruePositiveStd(ii) = mean(bigFatTreeTruePositiveTemp);
    bigFatTreeRelErrorsStd(ii) = mean(bigFatTreeRelErrorsTemp);
    bigFatTreeErrorsStd(ii) = mean(bigFatTreeErrorsTemp);

end
treeSeqFig = figure('Units', 'Inches', 'Position', [0, 0, 8, 10]);
ax1 = subplot(2,1,1);
hold on
h1Left = plot(distancesToPlot,treeRelErrorsMean);
%h1Left = plot(distancesToPlot,treeRelErrorsMean+treeRelErrorsStd./sqrt(nReps)'.*tinv(0.975,nReps-1)','k--');
%h1Left = plot(distancesToPlot,treeRelErrorsMean-treeRelErrorsStd./sqrt(nReps)'.*tinv(0.975,nReps-1)','k--');
hold on
yyaxis right
h1Right = plot(distancesToPlot,treeTruePositiveMean);
deltaD = (dMax-dMin)/nBins;
xlim([dMin-deltaD/2 dMax+deltaD/2]);
xl = xlim;
ax1.YAxis(2).Color = [0 0 0];
legend('Relative Error','Parasite Misclassification Rate'...
            ,'Location','NorthWest')
title('(a) Classifier Performance')
yyaxis left
ylabel('Error Relative to random model')
yyaxis right
ylabel('Parasite Misclassification Rate')
        
        
ax2 = subplot(2,1,2);
yyaxis left
hold on
spy(predsUsedBin','k.');
h2Grid = plot([0;82],repmat(1:12,2,1)','-','Color',[0.8 0.8 0.8]);
xlim(xl*200);
hold off
title('(b) Predictors Used')
yyaxis right
plot(distancesToPlot*100+1,meanWebSizes);
%
set(ax2.YAxis(1)...
                   ,'TickValues',1:12 ...
                   ,'TickLabels',fancyNames(1:12)...
                   ,'Color',[0 0 0]...
                   );
set(ax2.YAxis(2)...
    ,'Limits',[floor(min(meanWebSizes)),ceil(max(meanWebSizes))]...
    ,'Color',[0 0 0]...
    );
ylabel('Mean Web Size')
yyaxis left
hold on
spy(predsUsedBin','k.');
hold off

ax2.PlotBoxAspectRatioMode = 'auto';
ax2.XAxis.TickValues = ax1.XAxis.TickValues*100+1;
ax2.XAxis.TickLabels = ax1.XAxis.TickLabels;
ax2.YAxis(2).MinorTick='on';

%}
xlabel('Minimum Clustering Distance')
arrayfun(@(x) set(x,'FontName',figureFont),treeSeqFig.Children)

%[left bottom width height]
ax1Pos = ax1.Position;
ax1Top = ax1Pos(2) + ax1Pos(4);
ax2Pos = ax2.Position;
ax2Shrink = 0.2*ax2Pos(4);
ax2Pos(4) = 0.8*ax2Pos(4);
ax1Pos(2) = ax2Pos(2) + ax2Pos(4) + 0.05;
ax1Pos(4) = ax1Top - ax1Pos(2) + 0.03;
ax2.Position = ax2Pos;
ax1.Position = ax1Pos;
print(sprintf('../figures/ParasiteAcc%sLinkage.png',linkageType),'-dpng','-r0')

save(sprintf('../../Chapter3/code/selectors%sLinkge.mat',linkageType),'selectors');
%
%save(sprintf('predFrequency%sLinkge.mat',linkageType),'predFrequency');

%}
%
end
function [] = plotABar(m,s,ref)


hold on
h1 = bar(m);
q = refline(0,ref);
q.LineStyle = '--';
q.Color = [.7 .7 .7];
e1 = errorbar(1:numel(m),m,s/sqrt(100)*1.96);
e1.LineStyle = 'none';
e1.Color = 'k';
h1.FaceAlpha = 0.3;
hold off

end