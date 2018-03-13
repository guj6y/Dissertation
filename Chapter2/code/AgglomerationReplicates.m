
webs = 1:12;
nWebs = numel(webs);
props = cell(nWebs,1);

meanVarStructure = struct('sum',[],...
                          'sq',[],...
                          'mean',[],...
                          'var',[]);

propStructure = struct('global',meanVarStructure,...
                       'local',meanVarStructure);
                   
linkageStructure = struct('max',propStructure,...
                          'mean',propStructure,...
                          'min',propStructure);

[props{:}] = deal(linkageStructure);
nReps = 100;
for ii = webs
    %This replicates the agglomeration procedure; account for any
    %randomness in the agglomeration method.
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
    end
    
    
     localZeros = zeros(S,2,4);
     globalZeros = zeros(S,10);
    
    [props{ii}.max.global.sum,...
     props{ii}.max.global.sq,...
     props{ii}.mean.global.sum,...
     props{ii}.mean.global.sq,...
     props{ii}.min.global.sum,...
     props{ii}.min.global.sq] = deal(globalZeros);
    
     [props{ii}.max.local.sum,...
     props{ii}.max.local.sq,...
     props{ii}.mean.local.sum,...
     props{ii}.mean.local.sq,...
     props{ii}.min.local.sum,...
     props{ii}.min.local.sq] = deal(localZeros);
    
    for jj = 1:nReps
        if ii >= 7
            
            [res, con] = NicheModel_nk(S,C);
            consumers = unique(con);
            nCon = numel(consumers);
            para = false(S,1);
            nPar = round(fPar*nCon);
            paraIds = randsample(consumers,nPar);
            para(paraIds) = true;
        end
        [propMax, propMean, propMin] = Agglomerate(res,con,para);
        props{ii}.max.global.sum= props{ii}.max.global.sum+ propMax.global;
        props{ii}.max.global.sq= props{ii}.max.global.sq+ propMax.global.^2;
        
        props{ii}.max.local.sum = props{ii}.max.local.sum+propMax.local;
        props{ii}.max.local.sq = props{ii}.max.local.sq+propMax.local.^2;
        
        props{ii}.mean.global.sum= props{ii}.mean.global.sum+ propMean.global;
        props{ii}.mean.global.sq= props{ii}.mean.global.sq+ propMean.global.^2;
        
        props{ii}.mean.local.sum = props{ii}.mean.local.sum+propMean.local;
        props{ii}.mean.local.sq = props{ii}.mean.local.sq+propMean.local.^2;
        
        props{ii}.min.global.sum= props{ii}.min.global.sum+ propMin.global;
        props{ii}.min.global.sq= props{ii}.min.global.sq+ propMin.global.^2;
        
        props{ii}.min.local.sum = props{ii}.min.local.sum+propMin.local;
        props{ii}.min.local.sq = props{ii}.min.local.sq+propMin.local.^2;
    end
    
    props{ii}.max.global.mean = props{ii}.max.global.sum/nReps;
    props{ii}.max.global.var = props{ii}.max.global.sq/(nReps-1) - (nReps/(nReps - 1))*props{ii}.max.global.mean.^2;
    
    props{ii}.max.local.mean = props{ii}.max.local.sum/nReps;
    props{ii}.max.local.var = props{ii}.max.local.sq/(nReps-1) - (nReps/(nReps - 1))*props{ii}.max.local.mean.^2;
    
    props{ii}.mean.global.mean = props{ii}.mean.global.sum/nReps;
    props{ii}.mean.global.var = props{ii}.mean.global.sq/(nReps-1) - (nReps/(nReps - 1))*props{ii}.mean.global.mean.^2;
    
    props{ii}.mean.local.mean = props{ii}.mean.local.sum/nReps;
    props{ii}.mean.local.var = props{ii}.mean.local.sq/(nReps-1) - (nReps/(nReps - 1))*props{ii}.mean.local.mean.^2;
    
    props{ii}.min.global.mean = props{ii}.min.global.sum/nReps;
    props{ii}.min.global.var = props{ii}.min.global.sq/(nReps-1) - (nReps/(nReps - 1))*props{ii}.min.global.mean.^2;
    
    props{ii}.min.local.mean = props{ii}.min.local.sum/nReps;
    props{ii}.min.local.var = props{ii}.min.local.sq/(nReps-1) - (nReps/(nReps - 1))*props{ii}.min.local.mean.^2;

end

save('AgglomerationProps.mat','props')