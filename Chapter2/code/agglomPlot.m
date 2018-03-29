load('data/Processed/webGeneration.mat')
%load('data/Processed/webGenerationEctoInvert.mat')
webs = 1:6;
nWebs = numel(webs);



localNames = {'web'           ...1
                ,'vul'        ...2
                ,'gen'        ...3
                ,'patl'       ...4
                ,'cc'         ...5
                ,'pr'         ...6
                ,'btwn'       ...7
                ,'ecoBtwn'    ...8
                ,'meanVulPrey'...9
                ,'meanGenPred'...10
                ,'para'       ...11
                ,'free'       ...12
                ,'S'          ...13
                ,'C'          ...14
                ,'size'       ...14
                ,'ccCyc'      ...15
                ,'ccMid'      ...16
                ,'ccIn'       ...17
                ,'ccOut'      ...18
                ,'counter'    ...19
                };


globalNames = { 'web'           ...1
              ,'S'              ...2
              ,'Cf'             ...3
              ,'Cff'            ...4
              ,'Cpp'            ...5
              ,'Cfp'            ...6
              ,'Cpf'            ...7
              ,'fPar'           ...8
              ,'numLargest'     ...9
              ,'n'              ...10
              ,'top'            ...11
              ,'int'            ...12
              ,'bas'            ...13
              ,'herb'           ...14
              ,'omn'            ...15
              ,'C'              ...16
              ,'d_J'            ...17
              ,'fCarn'          ...18
              ,'Cp'             ...19
              ,'Sf'             ...20
              ,'counter'        ...
              };
localMeanNames = {
                 'web'...                     1
                ,'vulFree'...                 2
                ,'vulPara'...                 3
                ,'genFree'...                 4
                ,'genPara'...                 5
                ,'patlFree'...                6
                ,'patlPara'...                7
                ,'ccFree'...                  8
                ,'ccPara'...                  9
                ,'prFree'...                  10
                ,'prPara'...                  11
                ,'btwnsFree'...               12
                ,'btwnsPara'...               13
                ,'ecoBtwnsFree'...            14
                ,'ecoBtwnsPara'...            15
                ,'meanVulPreyFree'...         16
                ,'meanVulPreyPara'...         17
                ,'meanGenPredFree'...         18
                ,'meanGenPredPara'...         19
                ,'S'...                       20
                ,'C'...                       21
                ,'minDistance'...             22
                ,'ccCycFree'        ...       23
                ,'ccCycPara'        ...       24
                ,'ccMidFree'        ...       25
                ,'ccMidPara'        ...       26
                ,'ccInFree'        ...        27
                ,'ccInPara'        ...        28
                ,'ccOutFree'        ...       29
                ,'ccOutPara'        ...       30
                ,'counter'                 ...31
            };
for linkageTypes = {'Max','Min','Mean'}
    %linkageType = 'Max';
    linkageType = linkageTypes{:};
    if strcmp(linkageType,'Max')
        fLinkage = 1e-10;
    elseif strcmp(linkageType,'Mean')
        fLinkage = 0.5;
    elseif strcmp(linkageType,'Min')
        fLinkage = 1;
    else
        %default to maximum linkage.
        linkageType = 'Max';
        fLinkage = 1e-10;
    end
    
    nReps = 1;
    propsLocal = cell(nWebs,1);
    propsGlobal = cell(nWebs,1);
    propsLocalMeans = cell(nWebs,1);
    fullProps = cell(nWebs,1);
    reducedProps = cell(nWebs,1);
    for ii = webs
        ii
        %This replicates the agglomeration procedure; account for any
        %randomness in the agglomeration method (there is very little for the
        %empirical webs (could get a precise number...) so don't worry about 
        %that for now.
        if ii <7
            res = linkListCell{ii}(:,1);
            con = linkListCell{ii}(:,2);
            para = propertiesCell{ii,2};
            S = numel(para);
        elseif ii>=7
            res = linkListCell{ii-6}(:,1);
            con = linkListCell{ii-6}(:,2);
            S = max([res;con]);
            L = numel(res);
            C = L/S^2;
            nPar = sum(propertiesCell{ii-6,2});
            nCon = S - numel(setdiff(res,con));
            fPar = nPar/nCon;
            %ALWAYS need 100 repetitions of the niche models, to get an average
            %response.
            nReps = 100;
        end


        for jj = 1:nReps
            if ii >= 7

                [res, con] = NicheModel_nk(S,C);
                consumers = unique(con);
                nCon = numel(consumers);
                para = false(S,1);
                nPar = round(fPar*nCon);
                paraIds = randsample(consumers,nPar);
                para(paraIds) = true;

                %Do something about the replications.
            else
            [properties,fullProps{ii}, reducedProps{ii}] = agglom(res,con,para,fLinkage);
            propsLocal{ii} = [repmat(ii,length(properties.local),1),properties.local];
            propsGlobal{ii} = [repmat(ii,length(properties.global),1), properties.global];
            propsLocalMeans{ii} = [repmat(ii,length(properties.localMeans),1), properties.localMeans];
            end

        end
    end


    localCol = containers.Map(localNames, 1:numel(localNames));
    globalCol = containers.Map(globalNames, 1:numel(globalNames));
    localMeansCol = containers.Map(localMeanNames, 1:numel(localMeanNames));
    propsLocal = cell2mat(propsLocal);
    propsGlobal = cell2mat(propsGlobal);
    propsLocalMeans = cell2mat(propsLocalMeans);
    save(sprintf('AgglomerationProps%sLinkage.mat',linkageType)...
        ...'AgglomerationPropsEctoInvert.mat'
                ,'propsLocal'...
                ,'propsGlobal'...
                ,'propsLocalMeans'...
                ,'localCol'...
                ,'globalCol'...
                ,'localMeansCol'...
                ,'fullProps'...
                ,'reducedProps'...
                )
end