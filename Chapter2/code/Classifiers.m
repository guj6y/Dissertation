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

load('data/Processed/webGeneration.mat')
T = load('AgglomerationPropsMajority.mat','propsLocal','localCol','propsGlobal','globalCol');

propsLocal = T.propsLocal;
localCol = T.localCol;
propsGlobal = T.propsGlobal;
globalCol = T.globalCol;
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
nPred = numel(variableNames);


%propsThatCanBeNan = propsLocal(:,[10, 15,16,17,18]);
%obsThatAreNan = isnan(propsThatCanBeNan);
%propsThatCanBeNan(obsThatAreNan) = 0;
%propsLocal(:,[10, 15,16,17,18]) = propsThatCanBeNan;
%propsLocal = [propsLocal obsThatAreNan];

sMin = min(propsLocal(:,localCol('S')));
sMax = max(propsLocal(:,localCol('S')));

carnAndPara = (propsLocal(:,localCol('free'))>0.5)|(propsLocal(:,localCol('para'))>=0.5);
webs = 1:6;
selectThisLevelGlobal = arrayfun(@(x) find(propsGlobal(:,globalCol('web'))==x,1),webs);
Ss = propsGlobal(selectThisLevelGlobal,globalCol('S'));
clear propsGlobal globalCol
selectThisLevelLocal = arrayfun(@(x,y) (propsLocal(:,localCol('S'))==y)&(propsLocal(:,localCol('web'))==x),webs',Ss,'UniformOutput',false);

initialWebs = sum(cell2mat(selectThisLevelLocal'),2)>0;
initialWebsPropsLocal = propsLocal(initialWebs,:);
initCarn = (initialWebsPropsLocal(:,localCol('free'))>0.5);
initPara = (initialWebsPropsLocal(:,localCol('para'))>=0.5);
initCarnOrPara = initCarn|initPara;

X0 = propsLocal(carnAndPara,[2:4 6:10, 16:19]);
X0(isnan(X0)) = 0;

XForBig = initialWebsPropsLocal(initCarnOrPara,[2:4 6:10, 16:19]);
XForBig(isnan(XForBig)) = 0;
YForBig = initialWebsPropsLocal(initCarnOrPara,localCol('para'))>=0.5;
fParaBig = mean(YForBig);

nBig = numel(YForBig);
ogNodeList = 1:nBig;
nBigTest = 1000;
classificationPlot = figure;

impVars = [3 5 8 10 12];
XImp = XForBig(:,impVars);
varNames = variableNames(impVars);
impVars = [impVars 0];
nImpVars = numel(impVars);

relErrorMeans = zeros(nImpVars-1,2);
relErrorStds = zeros(nImpVars-1,2);
truePositiveMeans = zeros(nImpVars-1,2);
truePositiveStds = zeros(nImpVars-1,2);

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

                %{
                cutOffVals = 0:.01:1;
                bigGlm = fitglm(xTrain,yTrain...
                        ,'distribution','binomial'...
                        ,'varNames',varNames...
                        );
                YTrainHats = bigGlm.predict(xTrain)>=cutOffVals;
                cutOffBig = cutOffVals(find(sum(YTrainHats)<=fParaBig*numel(YForBig),1));
                %}

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
title('Relative Error, Removed Variable')
xlabel('Variable Removed')
ylabel('Error Relative to Random Model')


a2 = subplot(2,2,2);
plotABar(relErrorMeans(:,1),relErrorStds(:,1),fullRelErrorMean);
title('Relative Error, Isolated Variable')
xlabel('Variable Isolated')
ylabel('Error Relative to Random Model')


a3 = subplot(2,2,3);
plotABar(truePositiveMeans(:,2),truePositiveStds(:,2),fullTruePositiveMean);
title('True Positive, Removed Variable')
xlabel('Variable Removed')
ylabel('Fraction of Parasites Misclassifid')

a4 = subplot(2,2,4);
plotABar(truePositiveMeans(:,1),truePositiveStds(:,1),fullTruePositiveMean);
title('True Positive, Isolated Variable')
xlabel('Variable Isolated')
ylabel('Fraction of Parasites Misclassified')

arrayfun(@(x) set(x,'FontName','CMU Serif'),classFigure.Children);
%arrayfun(@(x) set(x,'Interpreter','LaTeX'),[a1.Children;a2.Children;a3.Children;a4.Children]);
arrayfun(@(x) set(x,'XTickLabels',varNames),classFigure.Children);
arrayfun(@(x) set(x,'XTickLabelRotation',20),classFigure.Children);
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 10]);
arrayfun(@(x) set(x,'XTick',1:5),classFigure.Children)
print('../figures/treeAnalysisc.png','-dpng','-r0')
%This is insanity. Hard enough to make one model; don't try to make an
%arbitrary number of them. They were all worthless; the big one is 
%worseless on the agglomerated webs, too; so why did I even bother with
%that?
%{

%Y01 = propsLocal(:,11);   %Continuous; for regression.
Y0 = propsLocal(:,localCol('para'))>=0.5; %categorical
sLevels = 30:109;
nLevels = numel(sLevels);



glms = cell(nLevels,1);
trees = cell(nLevels,1);

glmErrors = nan(nLevels,1);
treeErrors = nan(nLevels,1);
naddas  = nan(nLevels,1);

treeParasites = nan(nLevels,1);
glmParasites = nan(nLevels,1); 
treeTotals = nan(nLevels,1);
glmTotals = nan(nLevels,1);

bigGlmParasites = nan(nLevels,1);
bigTreeParasites = nan(nLevels,1);
bigGlmTotals = nan(nLevels,1);
bigTreeTotals = nan(nLevels,1);
bigGlmErrors = nan(nLevels,1);
bigTreeErrors = nan(nLevels,1);

for ii = 1:nLevels
    
    s= sLevels(ii)
    
    pickThese = (propsLocal(carnAndPara,localCol('S'))==s);
    
    X = X0(pickThese,:);
    %Y1 = Y01(pickThese,:); %continuous for regression?
    Y = Y0(pickThese,:);
    
    n = sum(pickThese);
    nTrain = round(0.5*n);
    nTest = n-nTrain;
    
    glmError = 0;
    treeError = 0;
    glmParasite = 0;
    treeParasite = 0;
    glmTotal = 0;
    treeTotal = 0;
    
    bigGlmError = 0;
    bigTreeError = 0;
    bigGlmParasite = 0;
    bigTreeParasite = 0;
    bigGlmTotal = 0;
    bigTreeTotal = 0;
    
    nadda = 0;
    nReps = 10;
    for z = 1:nReps
    trainThese = false(n,1);
    trainTheseIndices = randsample(n,nTrain);
    trainThese(trainTheseIndices) = true;
    testThese = ~trainThese;
    
    XTrain = X(trainThese,:);
    YTrain = Y(trainThese,:);
    
    XTest = X(testThese,:);
    YTest = Y(testThese,:);

    glm = fitglm(XTrain,YTrain...
        ,'distribution','binomial'...
        ,'varNames',variableNames);
    YTrainHats = glm.predict(XTrain)>=cutOffVals;
    fracPara = mean(YTrain);
    cutOff = cutOffVals(find(sum(YTrainHats)<=fracPara*numel(YTrain),1));
    
    YTestHat = (glm.predict(XTest)>=cutOff);
    glmError = (1 - sum(YTestHat == YTest)/nTest)/nReps + glmError;
    glmTotal = (sum(YTestHat)/nTest)/nReps + glmTotal;
    glmParasite = sum(YTest&YTestHat)/sum(YTest)/nReps + glmParasite;
    
    tree = fitctree(XTrain,YTrain...
        ,'PredictorNames',variableNames(1:nPred-1)...
        ,'MaxNumSplits',4 ...
        );
    YTestHat = tree.predict(XTest);
    treeError = (1 -sum(YTestHat == YTest)/nTest)/nReps + treeError;
    treeTotal = (sum(YTestHat)/nTest)/nReps + treeTotal;
    treeParasite = sum(YTest&YTestHat)/sum(YTest)/nReps + treeParasite;
    
    nadda = mean(YTest)/10 + nadda;
    
    YTestHat = treeBig.predict(XTest);
    bigTreeError = (1 -sum(YTestHat == YTest)/nTest)/nReps + bigTreeError;
    bigTreeTotal = (sum(YTestHat)/nTest)/nReps + bigTreeTotal;
    bigTreeParasite = sum(YTest&YTestHat)/sum(YTest)/nReps + bigTreeParasite;
    
    YTestHat = (glmBig.predict(XTest)<=cutOffBig);
    bigGlmError = (1 - sum(YTestHat == YTest)/nTest)/nReps + bigGlmError;
    bigGlmTotal = (sum(YTestHat)/nTest)/nReps + bigGlmTotal;
    bigGlmParasite = sum(YTest&YTestHat)/sum(YTest)/nReps + bigGlmParasite;
    
    end
    
    %lms{ii} = lm;
    glms{ii} = glm;
    trees{ii} = tree;
    
    glmParasites(ii) = glmParasite;
    treeParasites(ii) = treeParasite;
    glmTotals(ii) = glmTotal;
    treeTotals(ii) = treeTotal;
    glmErrors(ii) = glmError;
    treeErrors(ii) = treeError;
    
    bigGlmParasites(ii) = bigGlmParasite;
    bigTreeParasites(ii) = bigTreeParasite;
    bigGlmTotals(ii) = bigGlmTotal;
    bigTreeTotals(ii) = bigTreeTotal;
    bigGlmErrors(ii) = bigGlmError;
    bigTreeErrors(ii) = bigTreeError;
    
    naddas(ii) = nadda;
    
end
%}
%{
close all
figureFont = 'CMU Serif';
classFig = figure;
subplot(2,1,1);
naddas(naddas>0.5) = 1-naddas(naddas>0.5);
randomLoss = 2*(naddas.*(1-naddas));
h1 = plot(sLevels,[glmErrors./randomLoss treeErrors./randomLoss bigTreeErrors./randomLoss bigGlmErrors./randomLoss]);
q = refline(0,1);
q.LineStyle = '--';
q.Color = [0.8 0.8 0.8];
ylim([0 2]);
%If you want to look at how the big tree does; I don't think it's good
%because you have a lot of the same species at different points along the
%clustering sequence.
%glmTotals treeTotals bigGlmErrors bigTreeErrors bigGlmTotals bigTreeTotals])
legend('logit','tree','bigTree','bigLogit','Location','NorthEast')
%heheheh
%...   ,'Pred. Frac. Para logit','Pred. Frac. Para Tree','bigLogit','bigTree','pred Frac. bigLogit','pred Frac. bigTree')
title('(a) Loss for two classifiers')
xlabel('Web Size')
ylabel('Fraction of Random Classifier Error')


subplot(2,1,2);
h2 = plot(sLevels,[glmParasites, treeParasites naddas.*naddas bigTreeParasites bigGlmParasites]);
%Samesies about the big tree/glm stuff.
%, bigGlmParasites, bigTreeParasites])
legend('logit','tree','random','bigTree','bigGlm') %,'bigGlm','bigTree')
title('(b) Fraction of Parasite Identified Correctly')
xlabel('Web Size')
ylabel('Fraction of Parasites')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 8, 10])
q = refline(0,0.5);
q.LineStyle = '--';
q.Color = [0.8 0.8 0.8];
z = refline(0,mean(glmParasites,'omitnan'));
z.LineStyle = '--';
z.Color = [0 0 .5];
z = refline(0,mean(treeParasites,'omitnan'));
z.LineStyle = '--';
z.Color = [.5 0 0];
arrayfun(@(x) set(x,'FontSize',10),classFig.Children);
arrayfun(@(x) set(x,'FontName',figureFont),classFig.Children);
print('../figures/ParasiteAcc.png','-dpng','-r0')
%}
function [] = plotABar(m,s,ref)


hold on
h1 = bar(m);
q = refline(0,ref);
q.LineStyle = '--';
q.Color = [.7 .7 .7];
e1 = errorbar(1:5,m,s/sqrt(1000)*1.96);
e1.LineStyle = 'none';
e1.Color = 'k';
h1.FaceAlpha = 0.3;
hold off

end