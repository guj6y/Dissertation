if ~exist('simParams','var')
    load ../raw/simParams.mat
end

if ~exist('extcts','var')
    load ../raw/rawOutputs.mat
end

if ~exist('webData','var')
    load ../raw/webData.mat
end

if ~exist('fParAll0','var')
    load ../raw/metaSimData.mat
end

nSims = numel(extcts)
curWeb = 0;

starvedVec = cell(nSims,1);
fParVec = cell(nSims,1);
kFreeVec = cell(nSims,1);
kParaVec = cell(nSims,1);
concVec = cell(nSims,1);
refVec = cell(nSims,1);
noPreyVec = cell(nSims,1); 
nExtctVec = cell(nSims,1);
extctOrdVec = cell(nSims,1);
extctIdVec = cell(nSims,1);
webVec = cell(nSims,1);
genVec = cell(nSims,1);
paraVec = cell(nSims,1);
simNoVec = cell(nSims,1);


for ii = 1:nSims
    thisWeb = simParams{ii}.web;
    
    if curWeb ~= thisWeb
        
        curWeb = thisWeb;
        LL = webData{curWeb}.LL;
        res = LL(:,1);
        con = LL(:,2);
        A = full(sparse(res,con,1,40,40));

    end

    extct = extcts{ii};
    [~,nCols] = size(extct);
    nExtct = nCols -2;
    fracStarved = 0;
    fracNoPrey = 0;
    extct = extct(2:end,2:end-1);
    findExtctOrder = sum(extct==0,2);
    findExtctOrder = nExtct + 1 - findExtctOrder;

    totExtct = sum(findExtctOrder <= nExtct);
    
    starvedVec{ii} = zeros(totExtct,1);
    fParVec{ii} = zeros(totExtct,1);
    kFreeVec{ii} = zeros(totExtct,1);
    kParaVec{ii} = zeros(totExtct,1);
    concVec{ii} = zeros(totExtct,1);
    refVec{ii} = zeros(totExtct,1);
    noPreyVec{ii} = zeros(totExtct,1);
    nExtctVec{ii} = zeros(totExtct,1);
    extctOrdVec{ii} = zeros(totExtct,1);
    extctIdVec{ii} = zeros(totExtct,1);
    webVec{ii} = zeros(totExtct,1);
    genVec{ii} = zeros(totExtct,1);
    paraVec{ii} = zeros(totExtct,1);
    simNoVec{ii} = zeros(totExtct,1);

    extctSoFar = 0;
    for jj = 1:nExtct
        
        extctNow = findExtctOrder == jj;
         
        %Handling multiple extinction events. 
        nExtctNow = sum(extctNow);
         
        B = extct(:,jj);
        Bh = B.^1.2;
        dens = 0.5^1.2 + A'*Bh;
        yF = Bh(res)./(dens(con))*8;
        Fij = full(sparse(res,con,yF,40,40));
        extctConsumption = sum(Fij(:,extctNow));
        thisRange = (extctSoFar+1):(extctSoFar+nExtctNow);
        simNoVec{ii}(thisRange) = ii;
        starvedVec{ii}(thisRange) = sum(extctConsumption < 1);
        noPreyVec{ii}(thisRange) = sum(extctConsumption == 0);
        extctOrdVec{ii}(thisRange) = jj;
        extctIdVec{ii}(thisRange) = find(extctNow);
        nExtctVec{ii}(thisRange) = nExtct;
        fParVec{ii}(thisRange) = simParams{ii}.fPar;
        kFreeVec{ii}(thisRange) = simParams{ii}.kFree;
        kParaVec{ii}(thisRange) = simParams{ii}.kPara;
        concVec{ii}(thisRange) = simParams{ii}.fracPara;
        refVec{ii}(thisRange) = simParams{ii}.fracFree;
        webVec{ii}(thisRange) = simParams{ii}.web;
        genVec{ii}(thisRange) = sum(A(B>0,extctNow));
        paraVec{ii}(thisRange) = simParams{ii}.para(extctNow);
        
        extctSoFar = extctSoFar + nExtctNow; 
    end
end

starvedVec = cell2mat(starvedVec);
noPreyVec = cell2mat(noPreyVec);
extctOrdVec = cell2mat(extctOrdVec);
extctIdVec = cell2mat(extctIdVec);
nExtctVec = cell2mat(nExtctVec);
fParVec = cell2mat(fParVec);
kFreeVec = cell2mat(kFreeVec);
kParaVec = cell2mat(kParaVec);
concVec = cell2mat(concVec);
refVec = cell2mat(refVec);
genVec = cell2mat(genVec);
simNoVec = cell2mat(simNoVec);

