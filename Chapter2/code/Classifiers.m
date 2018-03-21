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
T = load('AgglomerationProps.mat','propsLocal','localCol');

propsLocal = T.propsLocal;
localCol = T.localCol;
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
                ,'size'       ...14
                ,'ccCyc'      ...15
                ,'ccMid'      ...16
                ,'ccIn'       ...17
                ,'ccOut'      ...18
                };
%}
variableNames = {'vul'        ...1
                ,'gen'        ...2
                ,'patl'       ...3
                ,'pr'         ...4
                ,'btwn'       ...5
                ,'ecoBtwn'    ...6
                ,'meanVulPrey'...7
                ,'meanGenPred'...8
                ,'ccCyc'      ...9
                ,'ccMid'      ...10
                ,'ccIn'       ...11
                ,'ccOut'      ...12
                ,'para'       ...13
                };
clear T

sMin = min(propsLocal(:,localCol('S')));
sMax = max(propsLocal(:,localCol('S')));


X0 = propsLocal(:,[2:4, 6:10, 15:18]);
%Y01 = propsLocal(:,11);   %Continuous; for regression.
Y0 = propsLocal(:,11)>0; %categorical
sLevels = 30:85;
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
    
pickForBig = propsLocal(:,localCol('S'))==100;
XForBig = X0(pickForBig,:);
YForBig = Y0(pickForBig,:);

 glmBig = fitglm(XForBig,YForBig...
        ,'distribution','binomial'...
        ,'varNames',variableNames);
   
 treeBig = fitctree(XForBig,YForBig...
        ,'PredictorNames',variableNames(1:12)...
        ,'MaxNumSplits',4 ...
        ,'MinleafSize',10);

for ii = 1:nLevels
    
    s= sLevels(ii);
    
    pickThese = (propsLocal(:,localCol('S'))==s);
    
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
    YTestHat = (glm.predict(XTest)>=.5);
    glmError = (1 - sum(YTestHat == YTest)/nTest)/nReps + glmError;
    glmTotal = (sum(YTestHat)/nTest)/nReps + glmTotal;
    glmParasite = sum(YTest&YTestHat)/sum(YTest)/nReps + glmParasite;
    
    tree = fitctree(XTrain,YTrain...
        ,'PredictorNames',variableNames(1:12)...
        ,'MaxNumSplits',4 ...
        ,'MinleafSize',10);
    YTestHat = tree.predict(XTest);
    treeError = (1 -sum(YTestHat == YTest)/nTest)/nReps + treeError;
    treeTotal = (sum(YTestHat)/nTest)/nReps + treeTotal;
    treeParasite = sum(YTest&YTestHat)/sum(YTest)/nReps + treeParasite;
    
    nadda = mean(YTest)/10 + nadda;
    
    YTestHat = treeBig.predict(XTest);
    bigTreeError = (1 -sum(YTestHat == YTest)/nTest)/nReps + bigTreeError;
    bigTreeTotal = (sum(YTestHat)/nTest)/nReps + bigTreeTotal;
    bigTreeParasite = sum(YTest&YTestHat)/sum(YTest)/nReps + bigTreeParasite;
    
    YTestHat = (glmBig.predict(XTest)>=.5);
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
close all
figureFont = 'CMU Serif';
lossFig = figure;
plot(sLevels,[glmErrors treeErrors naddas]);
%If you want to look at how the big tree does; I don't think it's good
%because you have a lot of the same species at different points along the
%clustering sequence.
%glmTotals treeTotals bigGlmErrors bigTreeErrors bigGlmTotals bigTreeTotals])
legend('logit','tree','frac. Para')
%heheheh
%...   ,'Pred. Frac. Para logit','Pred. Frac. Para Tree','bigLogit','bigTree','pred Frac. bigLogit','pred Frac. bigTree')
title('Loss for two classifiers')
xlabel('Web Size')
ylabel('Loss')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 17, 10])
arrayfun(@(x) set(x,'FontName',figureFont),lossFig.Children);
print('../figures/Loss.png','-dpng','-r0')
accFig = figure;
plot(sLevels,[glmParasites, treeParasites])
%Samesies about the big tree/glm stuff.
%, bigGlmParasites, bigTreeParasites])
legend('glm','tree') %,'bigGlm','bigTree')
title('Fractions Parasite Identified')
xlabel('Web Size')
ylabel('Fraction of Parasites')
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 17, 10])
arrayfun(@(x) set(x,'FontName',figureFont),accFig.Children);
print('../figures/ParasiteAcc.png','-dpng','-r0')