if ~exist('starvedVec','var')
    run checkExtinctions.m
end

if ~exist('fParAll0','var')
    load ../raw/metaSimData.mat
end

kFrees = [1 2 3];
kParas = [-3 -4 -5];
concs = [false true];
refs = [false true];

models = [1 2];
figures = [1 2 3];
cmap = jet(10);
for ii = figures
    figC = figure('Visible','off','Position',[0 0 800 600]);
    figP = figure('Visible','off','Position',[0 0 800 600]);
    
    count = 0;
    for jj = models
        for kk = kFrees 
            count = count + 1;
            maxN = -inf;
            
            for ll = 1:10
            selector = (kFreeVec == kk) & ...
                (concVec ==concs(jj))&(refVec==refs(jj)) &...
                (kParaVec == kParas(ii))&(fParVec == fParAll0(ll));
            plotThis = genVec(selector);
            figure(figC) 
            subplot(2,3,count);
            hold on
                [nThisModel,edgesThisModel] = histcounts(...
                    plotThis,...
                    'Normalization','cdf'...
                    ,'BinMethod','Integers');
                thisDist = line(edgesThisModel,[0,nThisModel]);
                set(thisDist,'Color',cmap(ll,:));
            hold off
            ylim([0 1])
            
            figure(figP) 
            subplot(2,3,count);
            hold on
                [nThisModel,edgesThisModel] = histcounts(...
                    plotThis,...
                    'Normalization','pdf'...
                    ,'BinMethod','Integers');
                thisDist = line(edgesThisModel,[0,nThisModel]);
                set(thisDist,'Color',cmap(ll,:));
                maxN = max(maxN,max(nThisModel));
            end
            hold off
            ylim([ 0 maxN])
        end
    end
    
    figure(figC);
    colormap(jet);
    caxis([0,0.5])
    colorbar('Position',[0.92 0.11 0.02 0.81],'Limits',[0,0.5]);
    print(figC...
        ,sprintf('../figures/perDistModel%u.jpg',ii)...
        ,'-djpeg','-r300')
    close(figC)
    
    figure(figP);
    colormap(jet);
    caxis([0,0.5])
    colorbar('Position',[0.92 0.11 0.02 0.81]);
    print(figP...
        ,sprintf('../figures/perDenModel%u.jpg',ii)...
        ,'-djpeg','-r300')
    close(figP)
end


