%This constructs the six web models I will look at.


%Correspond to .41, the most agglomeration before we start to really lose
%accuracy in the classification tree.
distancesPlotted = [0 0.05 0.1 0.2 0.4];

nCores = feature('numcores');
parpool('local',nCores);

for linkages = {'Min','Max'}

    linkageType = linkages{:};

    load(sprintf('../../Chapter2/code/selectors%sLinkage.mat',linkageType));

    %This is the average differences, std differences, etc. of all webs at each
    %agglomeration level. Also includes the diffs for each web, just in case
    %we want that, too.
    T = load(sprintf('groupedDistances%sLinkge.mat',linkageType));
    groupedDistanceMeans = T.groupedDistanceMeans;
    groupedDistanceStd = T.groupedDistanceStd;
    groupedDistanceN = T.groupedDistanceN;
    groupedDistanceTs = T.groupedDistanceTs;
    groupedDistancePs = T.groupedDistancePs;
    distancesToPlot = T.distancesToPlot;
    allMeans= T.plotByMinDistance;
    %We want the global properties only: need to get the Ss and Cs.
    T = load(sprintf('../../Chapter2/code/AgglomerationProps%sLinkage.mat',linkageType));

    propsGlobal = T.propsGlobal;
    globalCol = T.globalCol;
    adjMatrices = T.adjMatrices;
    propsLocal = T.propsLocal;
    localCol = T.localCol;

    T=load('../../Chapter2/code/data/Processed/webGeneration.mat','linkListCell');
    linkListCell = T.linkListCell;
    clear T

    adjMatrices = {...sparse(linkListCell{1}(:,1),linkListCell{1}(:,2),1,max(max(linkListCell{1})),max(max(linkListCell{1})))...
                  adjMatrices{1}{:}...
                  ...,sparse(linkListCell{2}(:,1),linkListCell{2}(:,2),1,max(max(linkListCell{2})),max(max(linkListCell{2})))...
                  ,adjMatrices{2}{:}...
                  ...,sparse(linkListCell{3}(:,1),linkListCell{3}(:,2),1,max(max(linkListCell{3})),max(max(linkListCell{3})))...
                  ,adjMatrices{3}{:}...
                  ...,sparse(linkListCell{4}(:,1),linkListCell{4}(:,2),1,max(max(linkListCell{4})),max(max(linkListCell{4})))...
                  ,adjMatrices{4}{:}...
                  ...,sparse(linkListCell{5}(:,1),linkListCell{5}(:,2),1,max(max(linkListCell{5})),max(max(linkListCell{5})))...
                  ,adjMatrices{5}{:}...
                  ...,sparse(linkListCell{6}(:,1),linkListCell{6}(:,2),1,max(max(linkListCell{6})),max(max(linkListCell{6})))...
                  ,adjMatrices{6}{:}...
                 }';

    
    webCodes = {'BSQ','CSM','EPB','FF','OH','STB'};
        localVariableNames = {'$v$'...          1
                             ,'$g$'...          2
                             ,'$T$'...         3    
                             ,'$Fr$'...           4
                             ,'$C_B$'...         5
                             ,'$C_{EB}$'...     6
                             ,'$v_r$'...  7
                             ,'$g_c$'...  8
                              ,'$\gamma^{cyc}$'...    14
                              ,'$\gamma^{mid}$'...    15
                              ,'$\gamma^{in}$'...     16
                              ,'$\gamma^{out}$'...    17
                             };
        globalVariableNames = {'TOP'...      1
                              ,'INT'...      2
                              ,'BAS'...      3
                              ,'HERB'...     4
                              ,'OMN'...      5
                              ,'CANN'...      6
                              ,'LOOP'...     7
                              ,'$\sigma_{Link}$'...   8
                              ,'$\sigma_g$'...    9
                              ,'$\sigma_v$'...    10
                              ,'TL'...       11
                              ,'maxSim'...   12
                              ,'$\mu_{PATH}$'...     13
                              ,'$\gamma^{cyc}$'...    14
                              ,'$\gamma^{mid}$'...    15
                              ,'$\gamma^{in}$'...     16
                              ,'$\gamma^{out}$'...    17
                              ,'$f_{carn}$'...    18
                           };
        modelNames = {'Random Consumer'...
                      ...,'Random Carnivore'...
                      ,'Tree Classified'...
                      ,'Inverse C'...
                      ,'Inverse Generalities'...
                      ,'Inverse Subwebs'};
    
    count = 0;
    
    for minDistance = distancesPlotted
        count = count+1;
        %first, need to get the empirical Data.
        %Agglomeration level:
        thisAggLevel = sum(distancesToPlot<=minDistance);

        sAndWeb = selectors(selectors(:,1) == thisAggLevel,[2 3]);
        getEm = sum((propsGlobal(:,globalCol('S'))==sAndWeb(:,1)' & propsGlobal(:,globalCol('web'))==sAndWeb(:,2)'),2)>0;
        getEmLocal = sum((propsLocal(:,localCol('S'))==sAndWeb(:,1)' & propsLocal(:,localCol('web'))==sAndWeb(:,2)'),2)>0;
        getEmLocalWebs = propsLocal(getEmLocal,localCol('web'));
        empPara = propsLocal(getEmLocal,localCol('para'));
        %Free is carn here...
        empCarn = propsLocal(getEmLocal,localCol('free'));
        useForTree = (empPara>=0.5)|(empCarn>=0.5);

        empParaCell = cellfun(@(x) empPara(getEmLocalWebs==x),{1,2,3,4,5,6},'UniformOutput',false);
        empWebs = adjMatrices(getEm);

        resCell = cell(6,1);
        conCell = cell(6,1);
        for ii = 1:6
            [resCell{ii},conCell{ii}] = ind2sub(size(empWebs{ii}),find(empWebs{ii}));
        end

        propertiesWebsToMatch = propsGlobal(getEm,:);


        %The plan is to calculate 100 webs with each proposed model for each web.
        %We then compare the average performance across all 6 webs.

        nWebsPerWeb = 1;

        webFunctions = {@nicheModelRandCon;
                        ...@nicheModelRandCarn;
                        @nicheModelRandTree;
                        @invNM1;
                        @invNM2;
                        @invNM3};
        nModels = numel(webFunctions);
        nWebs = 6;
        nLocalProps = 12;
        nGlobalProps = 18;

        carnParaDiffs = cell(nModels,nWebs);
        globalProps = cell(nModels,nWebs);
        treeProps = cell(nModels,nWebs,nWebsPerWeb);
        webFunctionOut = cell(nModels,nWebs,nWebsPerWeb);

        [carnParaDiffs{:}] = deal(zeros(nWebsPerWeb,nLocalProps));
        [globalProps{:}] = deal(zeros(nWebsPerWeb,nGlobalProps));

        %Get the empirical properties; after the next two lines,
        %empiricalproperties is 2x6; first row is global properties, second row is
        %local.
        empiricalProperties = cellfun(@(x,y,z) webPropsRaw(x,y,z>0),resCell,conCell,empParaCell','UniformOutput',false);
        empiricalProperties = reshape([empiricalProperties{:}],3,6);


        X = cell2mat(empiricalProperties(3,:)');
        X(isnan(X)) = 0;
        

        bigTree = fitctree(X(useForTree,:),empPara(useForTree)>=0.5 ...
                        ,'MaxNumSplits',4 ...
                        ,'PredictorNames',localVariableNames...
                            )


webProps = cell(nModels,nWebs);
        for ii = 1:nModels
           webFunction = webFunctions{ii};
           for jj = 1:nWebs
               %Just takes it all in; wanted this to be as standardized as
               %possible :)

               webFunctionIn = propertiesWebsToMatch(ii,:);
               tempGlobal = zeros(nWebsPerWeb,nGlobalProps);
               tempLocal = zeros(nWebsPerWeb,nLocalProps);
               parfor kk = 1:nWebsPerWeb
                   %WebFunctionOut is a cell array with the appropriate carn-para
                   %differences as well as the global properites.
                   webFunctionOut = webFunction(webFunctionIn,globalCol,bigTree);

                   tempGlobal(kk,:) = webFunctionOut{1};
                   tempLocal(kk,:) = webFunctionOut{2};
                   treeProps{ii,jj,kk} = webFunctionOut{3};
               end
               carnParaDiffs{ii,jj} = tempLocal;
               globalProps{ii,jj} = tempGlobal;
           end
        end
        save(sprintf('NicheTestResults%uDistance%sLinkageDistance%u',linkgeType,count)...
            ,'carnParaDiffs','globalProps','treeProps');
   end
end
clear carnParaDiffs globalProps treeProps
save('nicheModelTestsVariables.mat')
