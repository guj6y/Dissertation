%This is the average differences, std differences, etc. of all webs at each
%agglomeration level. Also includes the diffs for each web, just in case
%we want that, too.
%fileName = sprintf('NicheTestResults%sLinkage-Distance%u',linkageType,distID);
agglomData = sprintf('../../Chapter2/code/AgglomerationPropsMaxLinkage.mat');
load(agglomData)
linkageTypes = {'Max','Min'};

fracHitAllWebsEachModelG = zeros(5,4,2);
fracHitAllWebsEachModelF = zeros(5,4,2);
fracHitAllWebsEachModelP = zeros(5,4,2);

nPropsG = zeros(5,4,2);
nPropsF = zeros(5,4,2);
nPropsP = zeros(5,4,2);



fracNearMissedAllWebsEachModel = zeros(5,4,2);
fracMissedAllWebsEachModel = zeros(5,4,2);
nGlobalProps = 20;
nCommProps = 22;
nParaProps = nCommProps;
nFreeProps = nCommProps;
nModels = 4;
nWebs = 6;

distances = [0,0.05,0.1,0.2,0.4];



for linkageTypeID = 1:1
        
    linkageType = linkageTypes{linkageTypeID};
    T = load(sprintf('../../Chapter2/code/AgglomerationProps%sLinkage.mat',linkageType));
    propsGlobal = T.propsGlobal;
    globalCol = T.globalCol;
    reducedProps = T.reducedProps;

    for distanceID = 1
        nGlobalProps = 20;
        nCommProps = 22;
        minDistance = distances(distanceID);
        getEm = find(diff(propsGlobal(:,globalCol('d_J'))<=minDistance)==-1);
        propertiesWebsToMatch = propsGlobal(getEm,:);
        
        Ss = propertiesWebsToMatch(:,globalCol('S'));
        Cs = propertiesWebsToMatch(:,globalCol('C'));
        
        T = load(sprintf('NicheTestResults%sLinkage-Distance%u.mat',linkageType,distanceID));
        
        %We want the global properties only: need to get the Ss and Cs.
        S = load(sprintf('nicheModelTestsVariables.mat'));
        R = load('../../Chapter2/code/data/Processed/webGeneration.mat');

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
                           
            modelNames = {'Null Model'...
                          ...,'Random Carnivore'...
                          ...,'Tree Classified'...
                          ,'Inverse, C'...
                          ,'Inverse, $C_f$, $C_p$'...
                          ,'Inverse Niche Model '};

        %load data from NicheModelTests before doing this.

        MarkerCell = {'o';'+';'s';'*';'^';'x'};
        ColorCell = {'r';'b';'g';'c';'m';'k'};

        close all
        %globalFig = figure;
        %communityFig = figure;
        propsFig1 = figure('Visible','off');
        propsFig2 = figure('Visible','off');
        propsFig3 = figure('Visible','off');
        redundantGlobal = [3,5];
        redundantFree = [13,14];
        redundantPara = [13,14];

        
        useForPCg = true(1,nGlobalProps);
        useForPCg(redundantGlobal) = false;
        
        useForPCf = true(1,nCommProps);
        useForPCf(redundantFree) = false;

        %Need to find empirical properties >.<
        
        nDimsPredg = zeros(nModels,nWebs);
        nDimsPredf = zeros(nModels,nWebs);
        nDimsPredp = zeros(nModels,nWebs);
        
        %For each model
        nLinRedundG = zeros(nModels,nWebs);
        nLinRedundF = zeros(nModels,nWebs);
        nLinRedundP = zeros(nModels,nWebs);
        
        for ii = [1,4]
            
            useForPCp = true(1,nCommProps);
            useForPCp(redundantPara) = false;
            
            globalProps = T.globalProps(ii,:)';
            paraProps = T.paraProps(ii,:)';
            freeProps = T.freeProps(ii,:)';
            
            %For each web
            PCDataG = cell(size(globalProps));
            PCDataF = cell(size(paraProps));
            PCDataP = cell(size(freeProps));
            PCEmpProds = cell(3,6);
            empiricalProperties = S.empiricalProperties;
            
            for jj = 1:6
                %This is for PCA
                %nGlobalProps = 18;
                %nCommProps = 20;
                
                empThisWebG = empiricalProperties{1,jj}(:,useForPCg);
                empThisWebF = empiricalProperties{2,jj}(:,useForPCf);
                empThisWebP = empiricalProperties{3,jj}(:,useForPCp);

                thisWebGlobal = globalProps{jj}(:,useForPCg);
                thisWebFree = freeProps{jj}(:,useForPCf);
                thisWebPara = paraProps{jj}(:,useForPCp);
                
                VIFg = diag(inv(corr(thisWebGlobal)));
                VIFf = diag(inv(corr(thisWebFree)));
                VIFp = diag(inv(corr(thisWebPara)));
                
                goodG = VIFg<10;
                goodF = VIFf<10;
                goodP = VIFp<10;
                if (distanceID == 5)&& ii==4
                    fprintf('checkmeout\n')
                end
                nLinRedundG(ii,jj) = sum(~goodG);
                nLinRedundF(ii,jj) = sum(~goodF);
                nLinRedundP(ii,jj) = sum(~goodP);

                %{
                thisWebGlobal = thisWebGlobal(:,goodG);
                thisWebFree = thisWebFree(:,goodF);
                thisWebPara = thisWebPara(:,goodP);
                
                empThisWebG = empThisWebG(goodG);
                empThisWebF = empThisWebF(goodF);
                empThisWebP = empThisWebP(goodP);
                %}

                meansGlobal = mean(thisWebGlobal);
                meansFree = mean(thisWebFree);
                meansPara = mean(thisWebPara);

                stdGlobal = std(thisWebGlobal);
                stdFree = std(thisWebFree);
                stdPara = std(thisWebPara);

                okayG = stdGlobal > 0;
                okayF = stdFree > 0;
                okayP = stdPara > 0;
                
                empThisWebG = (empThisWebG(okayG) - meansGlobal(okayG))./stdGlobal(okayG);
                empThisWebF = (empThisWebF(okayF) - meansFree(okayF))./stdFree(okayF);
                empThisWebP = (empThisWebP(okayP) - meansPara(okayP))./stdPara(okayP);

                thisWebGlobal = (thisWebGlobal(:,okayG) - mean(thisWebGlobal(:,okayG)))./stdGlobal(:,okayG);
                thisWebPara = (thisWebPara(:,okayP) - mean(thisWebPara(:,okayP)))./stdPara(:,okayP);
                thisWebFree = (thisWebFree(:,okayF) - mean(thisWebFree(:,okayF)))./stdFree(:,okayF);

                [Ug, Sg, Vg] = svd(thisWebGlobal);
                [Up, Sp, Vp] = svd(thisWebPara);
                [Uf, Sf, Vf] = svd(thisWebFree);

                varEachDimg = cumsum(diag(Sg).^2)./sum(diag(Sg).^2);
                varEachDimf = cumsum(diag(Sf).^2)./sum(diag(Sf).^2);
                varEachDimp = cumsum(diag(Sp).^2)./sum(diag(Sp).^2);

                PCg = Ug*Sg;
                PCf = Uf*Sf;
                PCp = Up*Sp;
                
                PCDataG{jj} = PCg;
                PCDataF{jj} = PCf;
                PCDataP{jj} = PCp;
                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                empPCg = empThisWebG*Vg;
                empWebNormG = norm(empPCg);
                empPCgZ = empPCg./std(PCg);
                
                Pg = normcdf(abs(empPCgZ),'upper')*2;
                [PSortedg, PIdxg] = sort(Pg,'descend');
                PProdg = cumprod(PSortedg);
                dimPCg = numel(empPCgZ);
                sphereToCubeG = (pi^(dimPCg/2)/gamma(dimPCg/2+1)*empWebNormG^dimPCg)/prod(abs(2*empPCg));
                nDimsPredg(ii,jj) = sum(PProdg>.05./(1:dimPCg));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                empPCf = empThisWebF*Vf;
                empWebNormF = norm(empPCf);
                empPCfZ = empPCf./std(PCf);

                Pf = normcdf(abs(empPCfZ),'upper')*2;
                [PSortedf, PIdxf] = sort(Pf,'descend');
                PProdf = cumprod(PSortedf);
                dimPCf = numel(empPCfZ);
                sphereToCubeF = (pi^(dimPCf/2)/gamma(dimPCf/2+1)*empWebNormF^dimPCf)/prod(abs(2*empPCf));
                nDimsPredf(ii,jj) = sum(PProdf>.05./(1:dimPCf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                empPCp = empThisWebP*Vp;
                empWebNormP = norm(empPCp);
                empPCpZ = empPCp./std(PCp);

                Pp = normcdf(abs(empPCpZ),'upper')*2;
                [PSortedp, PIdxg] = sort(Pp,'descend');
                PProdp = cumprod(PSortedp);
                dimPCp = numel(empPCpZ);
                sphereToCubeP = (pi^(dimPCp/2)/gamma(dimPCp/2+1)*empWebNormP^dimPCp)/prod(abs(2*empPCp));
                nDimsPredp(ii,jj) = sum(PProdp>.05./(1:dimPCp));

                PCEmpProds{1,jj} = empPCg;
                PCEmpProds{2,jj} = empPCf;
                PCEmpProds{3,jj} = empPCp;
            end
            
        %{
        empiricalProperties = PCEmpProds;
        globalProps = PCDataG;
        freeProps = PCDataF;
        paraProps = PCDataP;
        nGlobalProps = sum(useForPCg);
        nFreeProps = sum(useForPCf);
        nParaProps = sum(useForPCp);
        %}
        %
        thisWebGlobal = globalProps; 
        thisWebFree = freeProps;
        thisWebPara = paraProps;
        %}
        globalPropsModEmp = cellfun(@(x,y) [x;y],globalProps',empiricalProperties(1,:),...
            'UniformOutput',false);
        freePropsModEmp = cellfun(@(x,y) [x;y],freeProps',empiricalProperties(2,:),...
            'UniformOutput',false);
        paraPropsModEmp = cellfun(@(x,y) [x;y],paraProps',empiricalProperties(3,:),...
            'UniformOutput',false);
        
        skewnessGlobal = cellfun(@skewness,globalProps,'UniformOutput',false);
        skewnessFree = cellfun(@skewness,freeProps,'UniformOutput',false);
        skewnessPara = cellfun(@skewness,paraProps,'UniformOutput',false);
        
        rightSkewG = cell2mat(cellfun(@(x) x>0.5,skewnessGlobal,'UniformOutput',false))';
        rightSkewF = cell2mat(cellfun(@(x) x>0.5,skewnessFree,'UniformOutput',false))';
        rightSkewP = cell2mat(cellfun(@(x) x>0.5,skewnessPara,'UniformOutput',false))';
        
        leftSkewG = cell2mat(cellfun(@(x) x<-0.5,skewnessGlobal,'UniformOutput',false))';
        leftSkewF = cell2mat(cellfun(@(x) x<-0.5,skewnessFree,'UniformOutput',false))';
        leftSkewP = cell2mat(cellfun(@(x) x<-0.5,skewnessPara,'UniformOutput',false))';
        
        [globalPropsModEmpSorted, globalPropsModEmpIdx] = ...
            (cellfun(@sort,globalPropsModEmp,'UniformOutput',false));
        [freePropsModEmpSorted, freePropsModEmpIdx] = ...
            (cellfun(@sort,freePropsModEmp,'UniformOutput',false));
        [paraPropsModEmpSorted, paraPropsModEmpIdx] = ...
            (cellfun(@sort,paraPropsModEmp,'UniformOutput',false));
        
        prctilesGlobal = cell2mat(cellfun(@(y) (find(y==1001)-0.5)/1001-(0:nGlobalProps-1)'...
            ,globalPropsModEmpIdx,'UniformOutput',false));
        
        prctilesFree = cell2mat(cellfun(@(y) (find(y==1001)-0.5)/1001-(0:nFreeProps-1)'...
            ,freePropsModEmpIdx,'UniformOutput',false));
        
        prctilesPara = cell2mat(cellfun(@(y) (find(y==1001)-0.5)/1001-(0:nParaProps-1)'...
            ,paraPropsModEmpIdx,'UniformOutput',false));
        
        prctilesPlotG = zeros(size(prctilesGlobal));
        prctilesPlotF = zeros(size(prctilesFree));
        prctilesPlotP = zeros(size(prctilesPara));
        
        prctilesPlotG(rightSkewG) = -log10(1-prctilesGlobal(rightSkewG));
        prctilesPlotF(rightSkewF) = -log10(1-prctilesFree(rightSkewF));
        prctilesPlotP(rightSkewP) = -log10(1-prctilesPara(rightSkewP));
        
        prctilesPlotG(leftSkewG) = log10(prctilesGlobal(leftSkewG));
        prctilesPlotF(leftSkewF) = log10(prctilesFree(leftSkewF));
        prctilesPlotP(leftSkewP) = log10(prctilesPara(leftSkewP));
        
        prctilesPlotG(~(leftSkewG|rightSkewG)) = log10(prctilesGlobal(~(leftSkewG|rightSkewG))) - log10(1-prctilesGlobal(~(leftSkewG|rightSkewG)));
        prctilesPlotF(~(leftSkewF|rightSkewF)) = log10(prctilesFree(~(leftSkewF|rightSkewF))) - log10(1-prctilesFree(~(leftSkewF|rightSkewF)));
        prctilesPlotP(~(leftSkewP|rightSkewP)) = log10(prctilesPara(~(leftSkewP|rightSkewP))) - log10(1-prctilesPara(~(leftSkewP|rightSkewP)));
       
        prctMax = log10(1000.5/1001) - log10(.5/1001);
        prctMin = log10(.5/1001) - log10(1000.5/1001);
        
        prct975 = log10(.975) - log10(.025);
        prct025 = -prct975;
        
        prct995 = log10(.995) - log10(.005);
        prct005 = -prct995;

        prct95 = log10(.95) - log10(.05);
        prct05 = -prct95;
       
        prct9875 = log10(.9875) - log10(.0125);
        prct0125 = -prct9875;
        

        nPtsG = numel(prctilesPlotG);
        nPtsF = numel(prctilesPlotF);
        nPtsP = numel(prctilesPlotP);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        prctOverG = (sum(prctilesPlotG(rightSkewG)>prct95) + ...
                   sum(prctilesPlotG(~(rightSkewG|leftSkewG))>prct975))...
                    /nPtsG;
                  
        prctOverF = (sum(prctilesPlotF(rightSkewF)>prct95) + ...
                   sum(prctilesPlotF(~(rightSkewF|leftSkewF))>prct975))...
                    /nPtsF;
            
        prctOverP = (sum(prctilesPlotP(rightSkewP)>prct95) + ...
                   sum(prctilesPlotP(~(rightSkewP|leftSkewP))>prct975))...
                    /nPtsP;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        prctUnderG = (sum(prctilesPlotG(leftSkewG)<prct05) + ...
                    sum(prctilesPlotG(~(rightSkewG|leftSkewG))<prct025))...
                    /nPtsG;
                  
        prctUnderF = (sum(prctilesPlotF(leftSkewF)<prct05) + ...
                    sum(prctilesPlotF(~(rightSkewF|leftSkewF))<prct025))...
                    /nPtsF;
            
        prctUnderP = (sum(prctilesPlotP(leftSkewP)<prct05) + ...
                    sum(prctilesPlotP(~(rightSkewP|leftSkewP))<prct025))...
                    /nPtsP;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        prctInG = 1 - prctOverG - prctUnderG;

        prctInF = 1 - prctOverF - prctUnderF;

        prctInP = 1 - prctOverP - prctUnderP;
        
        set(groot,'defaulttextinterpreter','latex'...
                 ,'defaulttextfontname','CMU Serif'...
                 ,'defaultaxesfontname','CMU Serif'...
            );
        
        %
        if distanceID==1
        xPara = repmat((1:6)',1,nParaProps);
        xFree = repmat((1:6)',1,nFreeProps);
        xGlobal = repmat((1:6)',1,nGlobalProps);
         
        figure(propsFig1);
        ax1 = subplot(1,2,round(sqrt(ii)));
        hold on
        a = scatter(xGlobal(:),prctilesPlotG(:));
        title(modelNames{ii})
        u1 = refline(0,prct975);    
        l1 = refline(0,prct025);
        txt1u = text(ax1,8,prct9875,sprintf('%.1f\\%%',prctOverG*100));
        txt1m = text(ax1,8,0,sprintf('%.1f\\%%',prctInG*100));
        txt1l = text(ax1,8,prct0125,sprintf('%.1f\\%%',prctUnderG*100));
        hold off

        figure(propsFig2);
        ax2 = subplot(1,2,round(sqrt(ii)));
        hold on
        b = scatter(xFree(:),prctilesPlotF(:));
        title(modelNames{ii})
        u2 = refline(0,prct975);
        l2 = refline(0,prct025);
        txt2u = text(ax2,8,prct9875,sprintf('%.1f\\%%',prctOverF*100));
        txt2m = text(ax2,8,0,sprintf('%.1f\\%%',prctInF*100));
        txt2l = text(ax2,8,prct0125,sprintf('%.1f\\%%',prctUnderF*100));
        hold off

        figure(propsFig3);
        ax3 = subplot(1,2,round(sqrt(ii)));
        hold on 
        c = scatter(xPara(:),prctilesPlotP(:)');
        title(modelNames{ii})
        u3 = refline(0,prct975);
        l3 = refline(0,prct025);
        txt3u = text(ax3,8,prct9875,sprintf('%.1f\\%%',prctOverP*100));
        txt3m = text(ax3,8,0,sprintf('%.1f\\%%',prctInP*100));
        txt3l = text(ax3,8,prct0125,sprintf('%.1f\\%%',prctUnderP*100));
        hold off

        arrayfun(@(x) set(x,'HorizontalAlignment','right'...
                           ,'FontSize',20 ...
                           ,'FontWeight','bold'...
                           ,'FontName','CMU Serif Bold')...
           ,[txt1u,txt1m,txt1l,txt2u,txt2m,txt2l,txt3u,txt3m,txt3l]...
                );
        
        if ii == 1
            ax1.YLabel.String = 'Entire Web ME';
            ax2.YLabel.String = 'Free Livers ME';
            ax3.YLabel.String = 'Parasites ME';
        end
        
        arrayfun(@(x) set(x,'XTick',1:6 ...
                           ,'XTickLabel',webCodes...
                           ,'XLim',[0,8]...
                           ,'YLim',[prctMin,prctMax]...
                           ,'FontName','CMU Serif'...
                           ,'FontSize',16 ...
                         )...
                 ,[ax1;ax2;ax3]);

        arrayfun(@(x) set(x,'MarkerFaceColor',[0 0 0]...
                           ,'MarkerFaceALpha',0.15 ...
                           ,'MarkerEdgeColor','none'...
                         )...
                ,[a;b;c]);
                     
        %arrayfun(@(x) set(x,'visible','off'),[b;c]);
         
         
        df = ones(64,1);

        hold off
        
        arrayfun(@(x) set(x,'LineStyle','--','Color',[0.1 0.1 0.1]...
                            ,'XData',[0,7]),[u1;l1;u2;l2;u3;l3]);
        set(propsFig1,'Units','Inches','Position',[0 0 17 10])
        set(propsFig2,'Units','Inches','Position',[0 0 17 10])
        set(propsFig3,'Units','Inches','Position',[0 0 17 10])
        end   
        
        %}
         %       
        allProps = prctilesPlotG(:);

        fracHitAllWebsEachModelG(distanceID,ii,linkageTypeID) = round(prctInG*nPtsG);
        nPropsG(distanceID,ii,linkageTypeID) = nPtsG; 
        
        fracHitAllWebsEachModelF(distanceID,ii,linkageTypeID) = round(prctInF*nPtsF);
        nPropsF(distanceID,ii,linkageTypeID) = nPtsF; 
        
        allProps = prctilesPlotP(:);
        fracHitAllWebsEachModelP(distanceID,ii,linkageTypeID) = round(prctInP*nPtsP);

        nPropsP(distanceID,ii,linkageTypeID) = nPtsP;
        %}
        end
    if distanceID == 1
        %print('../figures/Properties-PCA.png','-dpng','-r300')
        figure(propsFig1)
        print('../../Defense/figures/INM-Props1.png','-dpng','-r300')
        figure(propsFig2)
        print('../../Defense/figures/INM-Props2.png','-dpng','-r300')
        figure(propsFig3)
        print('../../Defense/figures/INM-Props3.png','-dpng','-r300')
    end
        
    end
end
%

LrsG = log10(...
    binopdf(fracHitAllWebsEachModelG(:,2:end,1),nPropsG(:,2:end,1),.95)./...
    binopdf(fracHitAllWebsEachModelG(:,1,1),nPropsG(:,1,1),.95));

LrsF = log10(...
    binopdf(fracHitAllWebsEachModelF(:,2:end,1),nPropsF(:,2:end,1),.95)./...
    binopdf(fracHitAllWebsEachModelF(:,1,1),nPropsF(:,1,1),.95));

LrsP = log10(...
    binopdf(fracHitAllWebsEachModelP(:,2:end,1),nPropsP(:,2:end,1),.95)./...
    binopdf(fracHitAllWebsEachModelP(:,1,1),nPropsP(:,1,1),.95));

modelPerfFig = figure('Visible','off','Units','Inches','Position',[0 0 11 14]);


ax1 = subplot(3,2,1);
p1= plot(distances',fracHitAllWebsEachModelG(:,[2,3,4,1],1)/nWebs);
legend('INM 1','INM 2','INM 3','NM 0')

ylabel('Avg. Global Props. Per Web')
xlabel('Min. Clust. Dist.','FontName','CMU Serif')

ax2 = subplot(3,2,2);
p2 = plot(distances',LrsG);

ylabel('Global Props. LR')
xlabel('Min. Clust. Dist.','FontName','CMU Serif')

ax3 = subplot(3,2,3);
p3 = plot(distances',fracHitAllWebsEachModelF(:,[2,3,4,1],1)/nWebs);
ylabel('Avg. FL Props. Per Web')
xlabel('Min. Clust. Dist.','FontName','CMU Serif')

ax4 = subplot(3,2,4);
p4 = plot(distances',LrsF);
ylabel('FL Props. LR')
xlabel('Min. Clust. Dist.','FontName','CMU Serif')

ax5 = subplot(3,2,5);
p5 = plot(distances',fracHitAllWebsEachModelP(:,[2,3,4,1],1)/nWebs);
ylabel('Avg. Frac. Para. Props. Per Web')
xlabel('Min. Clust. Dist.','FontName','CMU Serif')

ax6 = subplot(3,2,6);
p6 = plot(distances',LrsP);
ylabel('Para. Props. LR')
xlabel('Min. Clust. Dist.','FontName','CMU Serif')

arrayfun(@(x) set(x.Title,'String',sprintf('Average Number Properties per Web in 95\\%% Bounds'))...
         ,[ax1,ax3,ax5]);


arrayfun(@(x) set(x.Title,'String','Likelihood Ratio of Each model')...
         ,[ax2,ax4,ax6]);

arrayfun(@(x) set(x.Legend,'Location','best'...
                          ,'FontName','CMU Serif'...
                          ,'Interpreter','LaTeX'...
                 )...
         ,[ax2,ax4,ax6]);


arrayfun(@(x) set(x,'FontName','CMU Serif'...
               ,'XLim',[distances(1), distances(end)]...
               ,'FontName','CMU Serif'...
               ,'XTick',distances...
               ,'XTickLabels',cellfun(...
        @(y) sprintf('%.2g',y) ,mat2cell(distances,1,repmat(1,size(distances)))...
                                      ,'UniformOutput',false)...
               ,'FontName','CMU Serif'...
                 )...
         ,[ax1,ax2,ax3,ax4,ax5,ax6]...
        );
print('../figures/Fracs-LRs-PCA.png','-dpng','-r300')



 %}
