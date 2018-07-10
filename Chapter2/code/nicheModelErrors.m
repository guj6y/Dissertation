%This is the average differences, std differences, etc. of all webs at each
    %agglomeration level. Also includes the diffs for each web, just in case
    %we want that, too.
    %fileName = sprintf('NicheTestResults%sLinkage-Distance%u',linkageType,distID);
    S = load('nicheModelTestsVariables.mat');
    T = load('NicheTestResultsMaxLinkage-Distance1.mat');
    %We want the global properties only: need to get the Ss and Cs.

        localMeansCol = T.localMeansCol;
    localMeans = T.propsLocalMeans;

    webCodes = {'BSQ','CSM','EPB','FF','OH','STB'};
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
                  
%load data from NicheModelTests before doing this.

MarkerCell = {'o';'+';'s';'*';'^';'x'};
ColorCell = {'r';'b';'g';'c';'m';'k'};

carnParaSorted = cellfun(@sort,T.carnParaDiffs,'UniformOutput',false);
globalPropsSorted = cellfun(@sort,T.globalProps,'UniformOutput',false);
empiricalProperties = S.empiricalProperties;
close all
globalFig = figure;
communityFig = figure;
propsFig = figure;
for ii = 1:nModels
    prctsCarnPara = cellfun(@(x,y) log10((sum(x>=y)+1)/(nWebsPerWeb+2))-log10(1-(sum(x>=y)+1)/(nWebsPerWeb+2)),empiricalProperties(2,:),carnParaSorted(ii,:),'UniformOutput',false);
    prctsGlobalProp = cellfun(@(x,y) log10((sum(x>=y)+1)/(nWebsPerWeb+2))-log10(1-(sum(x>=y)+1)/(nWebsPerWeb+2)),empiricalProperties(2,:),carnParaSorted(ii,:),'UniformOutput',false);
    
    prctsCarnPara = cell2mat(prctsCarnPara');
    prctsGlobalProp = cel2mat(prctsGlobalProp');
    
    
    figure(globalFig);
    ax = subplot(nModels,1,ii);
    ax.XTick = 1:nGlobalProps;
    ax.XTickLabel = globalVariableNames;
    hold on
    hGlob = plot(prctsGlobalProp);
    hGlobCell = num2cell(hGlob);
    cellfun(@(x,y,z) set(x,'LineStyle','none'...
                          ,'Marker',y...
                          ,'Color',z...
                          ),hGlobCell,MarkerCell,ColorCell);
    axis([0.5, nGlobalProps+0.5, -4, 4]);
    sd2p = refline(0,2);
    sd2m = refline(0,-2);
    arrayfun(@(x) set(x,'LineStyle','--','Color',[0.7 0.7 0.7]),[sd2p;sd2m]);
    plot(repmat(1.5:(nGlobalProps-0.5),2,1),repmat([-4; 4],1,nGlobalProps-1),'Color',[0.8 0.8 0.8])
    sd0 = refline(0,0);
    sd0.Color = [0.7 0.7 0.7];
    title(modelNames{ii})
    hold off
    
    figure(communityFig);
    ax = subplot(nModels,1,ii);
    ax.XTick = 1:nGlobalProps;
    ax.XTickLabel = localVariableNames;
    hold on
    hLoc = plot(carnParaWebZ');
    hLocCell = num2cell(hLoc);
    cellfun(@(x,y,z) set(x,'LineStyle','none'...
                          ,'Marker',y...
                          ,'Color',z...
                          ),hLocCell,MarkerCell,ColorCell);
    axis([0.5, nLocalProps+0.5, -4, 4]);
    sd2p = refline(0,2);
    sd2m = refline(0,-2);
    arrayfun(@(x) set(x,'LineStyle','--','Color',[0.7 0.7 0.7]),[sd2p;sd2m]);
    plot(repmat(1.5:(nLocalProps-0.5),2,1),repmat([-4; 4],1,nLocalProps-1),'Color',[0.8 0.8 0.8])
    sd0 = refline(0,0);
    sd0.Color = [0.7 0.7 0.7];
    title(modelNames{ii})
    hold off
    
    figure(propsFig);
    ax = subplot(2,3,ii);
    xlim([0,7]);
    ax.XTick = 1:6;
    ax.XTickLabel = webCodes;
    hold on
    a = plot(carnParaWebZ,'ro');
    b = plot(globalWebZ,'k+');
    sd5u = refline(0,5);
    sd5l = refline(0,-5);
    arrayfun(@(x) set(x,'LineStyle','--','Color',[0.7 0.7 0.7]),[sd5u;sd5l]);
    hold off
    title(modelNames{ii})

end
