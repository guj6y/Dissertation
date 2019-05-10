S = load('data/Processed/Generation.mat','propertiesCell','linkListCell');
linkListCell = S.linkListCell;
propertiesCell = S.propertiesCell;  
paraCell = propertiesCell(:,2);

degreeDistributionData = zeros(6,10);
Ss = zeros(6,1);
Cs = zeros(6,1);
propertyCodes = {'genSD','vulSD','linkSD','corr_gv'};
   
localPropsPlottedFree = cellfun(@(x) strcat(x,'Free'),propertyCodes,'UniformOutput',false);
localPropsPlottedPara = cellfun(@(x) strcat(x,'Para'),propertyCodes,'UniformOutput',false);

propertyNames = {'S','C'...
                ,localPropsPlottedFree{:}...
                ,localPropsPlottedPara{:}...
                };
propertyDict = containers.Map(propertyNames, 1:length(propertyNames));

fancyLocalNames = {'(a) s_g'...
                  ,'(b) s_v'...
                  ,'(c) s_l'...
                  ,'(d) r_{gv}'};
fractionIneligable = zeros(4,12);

webs = 1:6;

webCodes = {'BSQ','CSM','EPB','FF','OH','STB'};
webCodes = webCodes(webs);
n_paras = zeros(size(webs));
n_carns = zeros(size(webs));
for i = 1:6
  ll = linkListCell{i};
  S = max(ll(:));
  Ss(i) = S;
  C = size(ll,1)/S^2;
  Cs(i) = C;
  para = paraCell{i};
  A = sparse(ll(:,1),ll(:,2),1,S,S);
  gen0 = sum(A)';
  vul0 = sum(A,2);
  # Lesson: Have different names for normalized and raw generalities/ vulnerabilities
  gen = gen0/mean(gen0);
  vul = vul0/mean(vul0);
  basal = gen == 0;
  carn = ~(basal) & ~(para) & (sum(A(basal,:))'==0);
  genSDPara = std(gen(para));
  vulSDPara = std(vul(para));
  linkSDPara = std(gen(para)+vul(para));
  corr_gvPara = corr(gen(para),vul(para));
  genSDFree = std(gen(carn));
  vulSDFree = std(vul(carn));
  linkSDFree = std(gen(carn) + vul(carn));
  corr_gvFree = corr(gen(carn),vul(carn));
  degreeDistributionData(i,:) = [S,C...
                                 ,genSDFree, vulSDFree...
                                 ,linkSDFree, corr_gvFree...
                                 ,genSDPara, vulSDPara...
                                 ,linkSDPara, corr_gvPara];
  n_elig_cc_c_para = sum(vul0(para)>1);
  n_elig_cc_c_free = sum(vul0(carn)>1);
  n_elig_cc_r_para = sum(gen0(para)>1);
  n_elig_cc_r_free = sum(gen0(carn)>1);
  n_elig_cc_cr_para = sum(vul0(para) > 0);
  n_elig_cc_cr_free = sum(vul0(carn) > 0);
  numberEligible(:,(i-1)*2 + 1) = [n_elig_cc_c_para;
                                   n_elig_cc_r_para;
                                   n_elig_cc_cr_para;
                                   n_elig_cc_cr_para];
  numberEligible(:,(i-1)*2 + 2) = [n_elig_cc_c_free;
                                   n_elig_cc_r_free;
                                   n_elig_cc_cr_free;
                                   n_elig_cc_cr_free];
  n_paras(i) = sum(para);
  n_carns(i) = sum(carn);
  fprintf('%s\t%u\t%.3f\t%u\t%u\n',webCodes{i}, C, numel(carn), sum(para), sum(carn));                                   
endfor



 
figureFont =  'Times New Roman';
alpha = 0.05;
m = 4;
alphaSeq = alpha./(m+1-(1:m));

textColRejected = [0.9 0.2 0.2];
textColFailed = [0.1 0.1 0.3];
lineColRejected = [1 0.5 0.5];
lineColFailed = [0.8 0.8 0.8];


n_elig_cell = mat2cell(numberEligible,repmat(1,4,1),repmat(1,12,1));
cc_codes = {'consumer','resource','consumer-resource','resource-consumer'}';
data_cell = horzcat(cc_codes,n_elig_cell);
fid = fopen('eligible_cc_counts.csv','w');
fprintf(fid,',%s,',webCodes{:});
fprintf(fid,'\nClustering Type')
ns = [n_paras; n_carns];
ns = ns(:);
n_cell = mat2cell(ns,repmat(1,12,1));
header = repmat({'n_para (Sp=','n_carn (Sc='},1,6);
header = cellfun(@(x,y) sprintf('%s%u)',x,y),header,n_cell','UniformOutput',false);
fprintf(fid,',%s',header{:})
fprintf(fid,'\n')
fmt = sprintf('%%s%s\\n',repmat(',%.5f',1,12));
fprintf(fid,fmt,data_cell'{:});
fclose(fid);

AR = 6/5.5;
img_width = 12;
img_height = img_width/AR;
degreeDistFig = figure('units','inches'...
                       ,'position',[0,0,img_width,img_height]); 
                       
colormap(summer);
cm_ = get(degreeDistFig,'colormap');
set(degreeDistFig,'colormap',cm_(size(cm_,1):-1:1,:));
set(degreeDistFig, 'paperunits', 'inches'...
                 , 'papertype','<custom>'...
                 , 'papersize', [img_width, img_height]...
                 , 'paperposition', [0, 0, img_width, img_height]...
    )
% Setup geometry of figure:
n_row = 2;
n_col = 2;

nLocalPropsPlotted = n_row*n_col;

left_offset = 0.5;
bot_offset = 0.04;
h_fig_sep = .02;
v_fig_sep = 0.02;
v_rem = 1 - bot_offset - n_row*v_fig_sep;
h_rem = 1 - .11 - h_fig_sep*(n_col - 1) - left_offset;
fig_width = h_rem/n_col;
fig_heighth = v_rem/n_row;

cMapMin = min(Cs)*.99-.001;
cMapMax = max(Cs)*1.01+.001;
colMatrix = summer(100);    
colMatrix(size(colMatrix,1):-1:1,:);
cVals = linspace(cMapMin,cMapMax,100);


tDiffs = zeros(nLocalPropsPlotted,1);
pDiffs = zeros(nLocalPropsPlotted,1);
xMeans = zeros(nLocalPropsPlotted,1);
yMeans = zeros(nLocalPropsPlotted,1);

rC = zeros(nLocalPropsPlotted,1);
rS = zeros(nLocalPropsPlotted,1);

for jj = 1:nLocalPropsPlotted
    figure(degreeDistFig);
    xs = degreeDistributionData(:,propertyDict(localPropsPlottedFree{jj}));
    ys = degreeDistributionData(:,propertyDict(localPropsPlottedPara{jj}));
    markSizes = sqrt(Ss-80)*30;
    markColors = colMatrix(arrayfun(@(x)find(x<cVals,1),Cs),:);
    fig_col = mod(jj-1,n_col) + 1;
    fig_row = floor((jj-1)/n_col)+1;
    h_pos = left_offset + (fig_row-1)*(fig_width + h_fig_sep);
    v_pos = bot_offset + (fig_col-1)*(fig_heighth + v_fig_sep);
    subplot(n_row,n_col,jj,'position',[h_pos, v_pos, fig_width, fig_heighth]);
    h = scatter(xs,ys,markSizes,markColors,'filled');
    rC(jj) = corr(ys-xs,Cs);
    rS(jj) = corr(ys-xs,Ss);
    if jj > n_col*(n_row-1)
        xlabel('Carnivore Score')
    else
        xlabel('')
    end
    if mod(jj,n_col) == 1
        ylabel('Parasite Score')
    else
        ylabel('')
    end
    title(fancyLocalNames{jj}...
        ,'Interpreter','tex'...
        ,'FontSize',12)
    xl = xlim;
    yl = ylim;
    xMeans(jj) = mean(xs);
    yMeans(jj) = mean(ys);
    diffs = ys-xs;
    meanDiff = mean(diffs);
    stdDiff = std(diffs)/sqrt(numel(diffs));
    tDiff = meanDiff./stdDiff;
    tDiffs(jj) = tDiff;
    pDiff = 1-tcdf(abs(tDiff),numel(diffs)-1);
    pDiffs(jj) = pDiff;
    lb = min(xl(1),yl(1));
    ub = max(xl(2),yl(2));
    
    %rl = refline(1,0);
    %rl.Color = textColFailed;
    % uistack(h,'top')
    hold on

    rl = line([lb ub],[lb ub],'color',textColFailed);
    h = scatter(xs,ys,markSizes,markColors);
    axis([lb ub lb ub]);
end

[pDiffsSorted,Idx] = sort(pDiffs);
[~,IdxIdx] = sort(Idx);
rejected = false(nLocalPropsPlotted,1);
rejected(Idx) = pDiffsSorted<alphaSeq';

tC = rC*2./sqrt(1-rC.^2);
tS = rS*2./sqrt(1-rS.^2);
pC = (1-tcdf(abs(tC),4))*2;
pS = (1-tcdf(abs(tS),4))*2;

CS_tab = [rC tC pC rS tS pS];
data_cell = horzcat(fancyLocalNames',mat2cell(CS_tab,repmat([1],n_row*n_col,1),repmat([1],6,1)));
fid = fopen('DegreeInteraction_cs.csv','w')
fprintf(fid,',corr(C),t(C),P(C),corr(S),t(S),P(S)\n')
fprintf(fid,'%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n',data_cell'{:})
fclose(fid)
for jj = 1:nLocalPropsPlotted

subplot(n_row,n_col,jj);
    hold on
    xl = xlim();
    yl = ylim();
    
    if rejected(jj)
        textCol = textColRejected;
        plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
        plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
        meanLabel = sprintf('t=%.3g\nP=%.3g<%.3g',tDiffs(jj),pDiffs(jj),alphaSeq(IdxIdx(jj)));
    else
        textCol = textColFailed;
        plot(xl,repmat(mean(yMeans(jj)),1,2),'Color',textCol);
        plot(repmat(mean(xMeans(jj)),1,2),yl,'Color',textCol);
        meanLabel = sprintf('t=%.3g\nP=%.3g\\geq%.3g',tDiffs(jj),pDiffs(jj),alphaSeq(IdxIdx(jj)));
    end
    
    textLocOpts = [xl(1)+diff(xl)*.02,yl(1)+diff(yl)*.02;...LowerLeft
                   xl(1)+diff(xl)*.02,yl(2)-diff(yl)*.02;...UpperLeft
                   xl(2)-diff(xl)*.02,yl(2)-diff(yl)*.02;...UpperRight
                   xl(2)-diff(xl)*.02,yl(1)+diff(yl)*.02]; %LowerRight
               
    alignmentOpts = {'VerticalAlignment','bottom','HorizontalAlignment','left';
                    'VerticalAlignment','top','HorizontalAlignment','left';
                    'VerticalAlignment','top','HorizontalAlignment','right';
                    'VerticalAlignment','bottom','HorizontalAlignment','right'};
    distToCentroids = sqrt(sum((textLocOpts-[xMeans(jj),yMeans(jj)]).^2,2));
    thisLoc = distToCentroids==max(distToCentroids);
    textLoc = textLocOpts(thisLoc,:);
    thisAlign = alignmentOpts(thisLoc,:);
    text(textLoc(1),textLoc(2),meanLabel...
                    ,'Color',textCol...
                    ,'FontSize',12 ...
                    ,thisAlign{:}...
                    ,'Margin',1 ...
                    ,'backgroundColor','white'...
                    ,'FontName',figureFont...
                    ,'Interpreter','TeX'...
                    ,'Fontweight','bold'...
                     );
    hold off
    axis([xl yl]);
end
a = colorbar('Position', [0.92 0.11 0.02 .81]...
             ,'Limits',[cMapMin,cMapMax]...
             );

cMapMid = mean([cMapMin,cMapMax]);
if min(abs(cMapMid-Cs))<.002
    cMapMid = [];
end
newTicks = [cMapMin; cMapMax; cMapMid; Cs];
[newTicks, Idx] = sort(newTicks);
% a.Ticks = newTicks;
webCodes = webCodes';
if isempty(cMapMid)
    newTickLabels = [{num2str(cMapMin);
                      num2str(cMapMax)};
                      webCodes(:)];  
else

    newTickLabels = [{num2str(cMapMin);
                      num2str(cMapMax);
                      num2str(cMapMid)};
                      webCodes(:)];
end

% a.TickLabels = newTickLabels(Idx);
set(a, 'ytick', newTicks);
set(a, 'yticklabel', newTickLabels(Idx));
caxis([cMapMin,cMapMax])
title(a,'C')

arrayfun(@(x) set(x,'FontName',figureFont),get(degreeDistFig,'children'));
fname = sprintf('../figures/InitialDegreeDist.png');
print (fname) -r600 -tight
