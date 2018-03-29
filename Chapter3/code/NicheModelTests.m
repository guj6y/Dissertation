%This constructs the six web models I will look at.


%Correspond to .41, the most agglomeration before we start to really lose
%accuracy in the classification tree.
distancesPlotted = [0 0.05 0.1 0.2 0.4];
addpath(genpath('~/matlab_bgl'));
%nCores = feature('numcores');
%parpool('local',nCores);

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
    reducedProps = T.reducedProps;
    propsLocal = T.propsLocal;
    localCol = T.localCol;
    localMeansCol = T.localMeansCol;
    localMeans = T.propsLocalMeans;

    clear T

    adjMatrices = {reducedProps{1}.linkage{:}...
                  ,reducedProps{2}.linkage{:}...
                  ,reducedProps{3}.linkage{:}...
                  ,reducedProps{4}.linkage{:}...
                  ,reducedProps{5}.linkage{:}...
                  ,reducedProps{6}.linkage{:}...
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
        
        getEm = find(diff(propsGlobal(:,globalCol('d_J'))<=minDistance)==-1);
        webCounter = propsGlobal(getEm,[globalCol('web'),globalCol('counter')]);
        getEmLocal = sum((propsLocal(:,localCol('web'))==webCounter(:,1)' & propsLocal(:,localCol('counter'))==webCounter(:,2)'),2)>0;
        empPara = cellfun(@(x) x.para(:)',reducedProps,'UniformOutput',false)';
        empPara = [empPara{:}]';
        empPara = empPara(getEm);
        empCarn = cellfun(@(x) x.carn(:)',reducedProps,'UniformOutput',false)';
        empCarn = [empCarn{:}]';
        empCarn = empCarn(getEm);
        
        useForTree = cellfun(@(x,y) x|y,empPara,empCarn,'UniformOutput',false);

        empParaCell = empPara;
        empPara = cell2mat(empParaCell);
        empParaCell = empParaCell;
        empCarn = cell2mat(empCarn);
        useForTree = cell2mat(useForTree);
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
        empiricalProperties = cellfun(@(x,y,z) webPropsRaw(x,y,z>0),resCell,conCell,empParaCell,'UniformOutput',false);
        empiricalProperties = reshape([empiricalProperties{:}],3,6);


        X = cell2mat(empiricalProperties(3,:)');
        X(isnan(X)) = 0;
        

        bigTree = fitctree(X(useForTree,:),empPara(useForTree)>=0.5 ...
                        ,'MaxNumSplits',4 ...
                        ,'PredictorNames',localVariableNames...
                            );


webProps = cell(nModels,nWebs);
        for ii = 1:nModels
           webFunction = webFunctions{ii};
           for jj = 1:nWebs
               %Just takes it all in; wanted this to be as standardized as
               %possible :)

               webFunctionIn = propertiesWebsToMatch(ii,:);
               tempGlobal = zeros(nWebsPerWeb,nGlobalProps);
               tempLocal = zeros(nWebsPerWeb,nLocalProps);
               for kk = 1:nWebsPerWeb
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
        save(sprintf('NicheTestResults%uDistance%sLinkageDistance%u',linkageType,count)...
            ,'carnParaDiffs','globalProps','treeProps');
   end
end
clear carnParaDiffs globalProps treeProps
save('nicheModelTestsVariables.mat')
