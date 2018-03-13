T = load('AgglomerationProps.mat');
props = T.props;
webs = 1:6; %6 empirical webs + averages for nichewebs with same S,C as these webs.
nWebs = numel(webs);


%Local properties, x2 since we have free and para averages here.
nLocalProps = 18;
localNames = props{ii}.(linkage).mean.local.Properties.VariableNames;
%Global properies; this is more complicated since we have two 
nGlobalProps = 15;
globalNames = props{ii}.(linkage).mean.global.Properties.VariableNames;
%Linkage Criterion: How links are assigned between agglomerated nodes. Can
%be 'max', 'mean', or 'min'.
linkage = 'max';


close all
%globalFig0 = figure();
%localFig0 = figure();
localFigFirst = figure();
xmin = inf;
ymin = inf;
xmax = 0;
ymax = 0;
webs = [3 2 6 1 5 4];
sizes = zeros(6,1);
count = 0;
ys = zeros(nWebs,9);
xs = zeros(nWebs,9);
plotAtThisLevel = 75;
for ii = webs
    x = props{ii}.(linkage).mean.global.nNodes; 
    localFigAll = figure;
for jj = 1:2:nLocalProps
        
        
    
        figure(localFigAll);
        subplot(3,3,(jj+1)/2);
        hold on
            title(localNames{jj})
            plot(x,ones(size(x)),'k:');
        plot(x,props{ii}.(linkage).mean.local.(localNames{jj}),'b')
        plot(x,props{ii}.(linkage).mean.local.(localNames{jj+1}),'r')
        plot(x,props{ii+6}.(linkage).mean.local.(localNames{jj}),'k-')
        plot(x,props{ii+6}.(linkage).mean.local.(localNames{jj+1}),'k--')
        hold off
        
        figure(localFigFirst);
        subplot(3,3,(jj+1)/2);
        hold on
        title(localNames{jj})
        x0 = props{ii}.(linkage).mean.local.(localNames{jj})(x==plotAtThisLevel);
        y0 = props{ii}.(linkage).mean.local.(localNames{jj+1})(x==plotAtThisLevel); 
        plot(x0,y0,'.','MarkerSize',15,'MarkerSize',500*props{ii}.(linkage).mean.global.Ccon(x==plotAtThisLevel));
        xs(ii,(jj+1)/2) = x0;
        ys(ii,(jj+1)/2) = y0;
        xlabel('free')
        ylabel('para')
        xl= xlim;
        yl = ylim;
        xlim([0 xl(2)]);
        ylim([0 yl(2)]);
        hold off
        
end
    
    
    sortedSizes = sort(sizes);
    
    
    count = 0;
   
end
%
figure(localFigFirst);
diffs = ys-xs;
meanDiffs = mean(diffs);
stdDiffs = std(diffs);
tDiffs = meanDiffs./(stdDiffs/(sqrt(6)));
Ps = 2*tcdf(abs(tDiffs),5,'upper');
[PsSorted, Idx] = sort(Ps);
alpha = 0.05;
m = 9;
alphas = alpha./(m+1-(1:m));
rejected = PsSorted<alphas;
rejectedUnSorted = rejected(Idx);
for jj = 1:9
        
        subplot(3,3,jj);
        hold on
        xl = xlim();
        yl = ylim();
        refline(1,0);
        plot(xl,repmat(mean(ys(:,jj)),1,2),'k--');
        plot(repmat(mean(xs(:,jj)),1,2),yl,'k--');
        
        
        hold off
        axis([xl yl]);
end
    %}
for ii = 1:nWebs
    globalFigAll = figure();
    x = props{ii}.(linkage).mean.global.nNodes; 
count = 0;    
 for jj = [3,8,4,11,5,12,6,13,7,14]
        
        count = count+1;
        subplot(5,2,count)
        hold on
            title(globalNames{jj})
        plot(x,props{ii}.(linkage).mean.global.(globalNames{jj}),'b')
        plot(x,props{ii+6}.(linkage).mean.global.(globalNames{jj}),'k')
        hold off
 end
end
