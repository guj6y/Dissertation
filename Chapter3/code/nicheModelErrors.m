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

    for distanceID = 1:5
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
                           
            modelNames = {'Random Consumer'...
                          ...,'Random Carnivore'...
                          ...,'Tree Classified'...
                          ,'Inverse, C'...
                          ,'Inverse, $C_f$, $C_p$'...
                          ,'Inverse, all $C$''s'};

        %load data from NicheModelTests before doing this.

        MarkerCell = {'o';'+';'s';'*';'^';'x'};
        ColorCell = {'r';'b';'g';'c';'m';'k'};

        close all
        %globalFig = figure;
        %communityFig = figure;
        propsFig = figure;
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
        
        for ii = 1:4
            
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
        %{
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
        
        set(groot,'defaulttextinterpreter','latex')
        
        
        
        %{
        figure(propsFig);
        xPara = repmat((1:6)',1,nParaProps);
        xFree = repmat((1:6)',1,nFreeProps);
        xGlobal = repmat((1:6)',1,nGlobalProps);
        
        ax1 = subplot(3,4,ii);
        hold on
        a = scatter(xPara(:),prctilesPlotP(:));
        title(modelNames{ii})
        u1 = refline(0,prct975);    
        l1 = refline(0,prct025);
        hold off
        ax2 = subplot(3,4,ii+4);
        hold on
        b = scatter(xFree(:),prctilesPlotF(:));
        u2 = refline(0,prct975);
        l2 = refline(0,prct025);
        hold off
        ax3 = subplot(3,4,ii+8);
        hold on
        c = scatter(xGlobal(:),prctilesPlotG(:)');
        u3 = refline(0,prct975);
        l3 = refline(0,prct025);
        hold off
        if ii == 1
            ax1.YLabel.String = 'Entire Web ME';
            ax2.YLabel.String = 'Free Livers ME';
            ax3.YLabel.String = 'Parasites ME';
        end
        arrayfun(@(x) set(x,'XTick',1:6 ...
                           ,'XTickLabel',webCodes...
                           ,'XLim',[0,7]...
                           ,'YLim',[prctMin,prctMax]...
                           ,'FontName','CMU Serif'...
                         )...
                 ,[ax1;ax2;ax3]);

        arrayfun(@(x) set(x,'MarkerFaceColor',[0 0 0]...
                           ,'MarkerFaceALpha',0.08 ...
                           ,'MarkerEdgeColor','none'...
                         )...
                ,[a;b;c]);
                     
        %arrayfun(@(x) set(x,'visible','off'),[b;c]);
        
        
        df = ones(64,1);

        hold off
        
        arrayfun(@(x) set(x,'LineStyle','--','Color',[0.7 0.7 0.7]...
                            ,'XData',[0,7]),[u1;l1;u2;l2;u3;l3]);
        set(propsFig,'Units','Inches','Position',[0 0 17 10])
        
        
        %}
                
        allProps = prctilesPlotG(:);
        fracHitAllWebsEachModelG(distanceID,ii,linkageTypeID) = sum((allProps<=prct975)&(allProps>=prct025));
        
        allProps = prctilesPlotF(:);
        fracHitAllWebsEachModelF(distanceID,ii,linkageTypeID) = sum((allProps<=prct975)&(allProps>=prct025));
        
        allProps = prctilesPlotP(:);
        fracHitAllWebsEachModelP(distanceID,ii,linkageTypeID) = sum((allProps<=prct975)&(allProps>=prct025));
        
        end
        %print('../figures/Properties-Raw.png','-dpng','-r0')
    end
end




