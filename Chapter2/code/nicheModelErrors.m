%load data from NicheModelTests before doing this.

MarkerCell = {'o';'+';'s';'*';'^';'x'};
ColorCell = {'r';'b';'g';'c';'m';'k'};

close all
globalFig = figure;
communityFig = figure;
propsFig = figure;
for ii = 1:nModels
    prctsCarnPara = cellfun(@(x,y) log10((sum(x>=y)+1)/(nWebsPerWeb+2))-log10(1-(sum(x>=y)+1)/(nWebsPerWeb+2)),empiricalProperties(2,:),carnParaSorted(ii,:),'UniformOutput',false);
    prctsGlobalProp = cellfun(@(x,y) log10((sum(x>=y)+1)/(nWebsPerWeb+2))-log10(1-(sum(x>=y)+1)/(nWebsPerWeb+2)),empiricalProperties(2,:),carnParaSorted(ii,:),'UniformOutput',false);
    
    prctsCarnPara = cell2mat(prctsCarnPara');
    
    cellfun(@(x,y) log10((sum(x<=y)+1)/(nWebPerWeb+1))-log10(1-(sum(x<=y)+1)/(nWebPerWeb+1)),empiricalProperties(2,:),carnParaSorted(ii,:));
    globalWebZ = cellfun(@(x,y,z) (x - y)./z ,empiricalProperties(1,:),globalPropsMeans(ii,:),globalPropsStds(ii,:),'UniformOutput',false);
    globalWebZ = cell2mat(globalWebZ');
    
    figure(globalFig);
    ax = subplot(nModels,1,ii);
    ax.XTick = 1:nGlobalProps;
    ax.XTickLabel = globalVariableNames;
    hold on
    hGlob = plot(globalWebZ');
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
