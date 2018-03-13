function [TS,extcts,bs,rs] = integrateExperiment(simParams,params,S,nSims)

%This automatically makes the parallel pool use every core possible. on ocelote, their
%default is (for some reason?) 12.
nCores = feature('numcores');
parpool('local',nCores);

%initialize output arrays.
TS= zeros(S,1000,nSims);
extcts = cell(nSims,1);
bs = zeros(nSims,2);
rs = zeros(nSims,2);

parfor  (simNo = 1:nSims, nCores)
   
    p = params;
    simParam = simParams{simNo};
    
    %%% Define properties for this model
    kFree = simParam.kFree;
    kPara = simParam.kPara;
    modelCode = simParam.modelCode;
    
    LL = simParam.LL;
    
    res = LL(:,1);
    con = LL(:,2);
    
    patl = simParam.patl;
    basal = simParam.gr>0;
    
    
    %Assimilation efficiency is given in Brose et al as .45 and .85 for
    %herbivores and carnivores, respectively.  This would correspond to
    %the links we talkin bout.  Easy enough to figure:
    eij = zeros(size(res));
    eij(basal(res)) = p.eijBas;
    eij(~basal(res)) = p.eijCon;  %Future: add an eij for parasitic links
    wij = true(size(res));
    
    %Set the properties that we have; web-level properties, independent
    %of parasite identity.
    p.B0 = simParam.B0;
    
    p.res = res;
    p.con = con;
    p.eij = eij;
    p.wij = wij;
    p.basal = basal;
    p.r = simParam.gr;
    p.modelCode = modelCode;
    
    ZFree = 10^kFree;
    
    %%% Define the parameters that change with parasite identities.
    ax = zeros(S,1);
    y = zeros(S,1);
    M = zeros(S,1);
    para = simParam.para;
    free = ~para;
    
    %assimilation rates
    y(free) = p.yFree;
    y(para) = p.yPara;
    
    %scaling constant
    ax(free) = p.axFree;
    ax(para) = p.axPara;
    
    %bodymass
    M(free) = ZFree.^(patl(free)-1);
    M(para) = 10.^(kPara + kFree*(patl(para)-2));
    M(basal) = p.r(basal);
    %metabolic rate
    x = ax.*M.^(-0.25);
    x(basal) = 0;
    
    %Save parameters to parameter structre.
    p.x = x;
    p.para = para;
    p.M = M;
    p.yij = y(con);
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Integrate!
    sol = integrateParasiteExperiments(p);    
 
    %Save the extinction times and B0 snapshots. Get the intial conditions and the final conditions as well. Doesn't matter,r eally! 
    extcts{simNo} = [sol.extctTime;sol.B0s];
    endSol = sol.(sprintf('sol%u',sol.n));
    xEnd = endSol.x(end);
    xRange = linspace(xEnd,xEnd-999,1000);
    TS(:,:,simNo) = deval(endSol,xRange);
    %sol is too large to save. a single run can be as much as 77MB. Which
    %obviously doesn't scale well for 21000 runs.

    %Calculate abundance-body size curves here. Could shunt this off to a 
    %seperate file, but I think it is best to put it in here.
    abundance = mean(TS(:,:,simNo),2)./M;
    bsClass10 = round(log10(M));
    x10 = min(bsClass10):max(bsClass10);
    y10 = log10(sum(abundance.*(bsClass10==x10)));
    
    finiteY10 = isfinite(y10);
    x10 = x10(finiteY10);
    y10 = y10(finiteY10);
    
    r10 = corr(x10',y10','rows','complete');
    b10 = r10*std(y10,'omitnan')./std(x10);

    bsClass = round(log(M));
    x = min(bsClass):max(bsClass);
    y = log(sum(abundance.*(bsClass==x)));
    
    finiteY = isfinite(y);
    x = x(finiteY);
    y = y(finiteY);
    
    r = corr(x',y','rows','complete');
    b = r*std(y,'omitnan')./std(x);
    
    bs(simNo,:) = [b10,b];
    rs(simNo,:) = [r10^2,r^2];

end

end
