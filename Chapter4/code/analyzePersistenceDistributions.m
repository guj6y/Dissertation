%if ~exists('starvedVec','v')
%    run checkExtinctions.m
%end

if ~exist('persistences','var')
    load ../raw/out.mat
end

%Making plots of cdfs of persistences. This bit makes a 2x2 plot for each value
%of kPara. row: using model, col: kFree
%

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
            selector = (kFreeVec == kk) & ...
                (concVec ==jj)&(refVec==jj) &...
                (kParaVec == ii);
            persistencesThisModel = activities.free(:,:,kk,ii,jj,jj)...
                + activities.para(:,:,kk,ii,jj,jj);
            
            figure(figC) 
            subplot(2,3,count);
            hold on
            for ll = 1:10
                [nThisModel,edgesThisModel] = histcounts(...
                    persistencesThisModel(:,ll),...
                    'Normalization','cdf');
                thisDist = line(edgesThisModel,[0,nThisModel]);
                set(thisDist,'Color',cmap(ll,:));
            end
            hold off
            ylim([0 1])
            figure(figP) 
            subplot(2,3,count);
            hold on
            maxN = -inf;
            for ll = 1:10
                [nThisModel,edgesThisModel] = histcounts(...
                    persistencesThisModel(:,ll),...
                    'Normalization','pdf');
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



