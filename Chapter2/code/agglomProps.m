function [localProperties,globalProperties,localMeans] = agglomProps(A,para,carn,nodeSize,minDistance)

%Check for the connected component; don't worry about the web having basal
%species, etc. since we are not trying to develop a dynamical model for
%these; we would deal with that in the modeling phase anyway. We do want to
%get rid of disconnected elements, though. This should (hopefully) have the
%added benefit that we don't get as many warnings when calculating the
%trophic leve.
S0 = length(A);
[n,grps] = graphconncomp(sparse(A),'weak',true,'directed',true);

maxCompSize = 0;
compSizes = zeros(1,n);
for ii = 1:n
    sizeii = sum(grps==ii);
    compSizes(ii) = sizeii;
    if sizeii>maxCompSize
        maxCompSize=sizeii;
        largeComp = ii;
    end
end
numLargest = sum(compSizes==maxCompSize);

if numLargest>1
    if maxCompSize<=(S0/4)
        fprintf('Web has lost more than half of its nodes. Returning nans.\n')
        localProperties = nan(1,18);
        globalProperties = nan(1,15);
        return
    end
    comps = 1:n;
    fprintf('Tie in largest component... averaging them.\n')
    localProperties = zeros(1,9);
    globalProperties = zeros(1,15);
    %TODO: fix this with new structures.
    for ii = comps(compSizes==maxCompSize)
        B = A(grps==ii,grps==ii);
        para_ = para(grps==ii);
        [loc,glob] = agglomProps(B,para_);
        localProperties = localProperties+loc/numLargest;
        globalProperties = globalProperties+glob/numLargest;
    end
    globalProperties(end) = n;
    return
end
%I like winners. Winners get to stay.
winners = grps==largeComp;
A = A(winners,winners);
[res,con] = find(A);
para = para(winners);


%Now, we calculate properties using the largest connected component.
S = length(A);
L = sum(sum(A));

%In degree (generality)
gen0 = full(sum(A))';


%This defines parasite cutoffs; without this line we calculate weighted
%averages:
binPara = para>0.5;
basal = gen0==0;
binFree = ~binPara;
binFree(basal) = false;
binCons = (~basal);
%Identify carnivores Only!
binCarns = true(S,1);
binCarns(con(basal(res))) = false;

free = 1-para;
free(basal) = 0;

para = para.*nodeSize;
%Pick the appropriate other comparison.
%free = free*nodeSize;
free = carn.*nodeSize;

gen = gen0./(mean(gen0));
genFree = wmean(gen,free)/mean(gen(binCons));
genPara = wmean(gen,para)/mean(gen(binCons));

%Mean generality of predators
%Mean topological fraction of predators' diet
meanGenPred = zeros(S,1);
for ii = 1:S
    pred_ii = con(res==ii);
    genPred_ii = gen(pred_ii);
    meanGenPred(ii) = mean(genPred_ii);
    if numel(genPred_ii)==0
        meanGenPred(ii) = 0;
    end
end
meanGenPred = meanGenPred./mean(meanGenPred);
meanGenPredFree = wmean(meanGenPred,free)/mean(meanGenPred(binCons));
meanGenPredPara = wmean(meanGenPred,para)/mean(meanGenPred(binCons));

%out degree (vulnerability)
vul0 = full(sum(A,2));
vul = vul0./mean(vul0);
vulFree = wmean(vul,free)/mean(vul(binCons));
vulPara = wmean(vul,para)/mean(vul(binCons));


%Mean vulnerability of prey
%Mean topological fraction of preys' binConsumers
meanVulPrey = zeros(S,1);
for ii = 1:S
    prey_ii = res(con==ii);
    vulPrey_ii = vul(prey_ii);
    meanVulPrey(ii) = mean(vulPrey_ii);
    if numel(vulPrey_ii)==0
        meanVulPrey(ii) = 0;
    end
end

meanVulPrey = meanVulPrey/mean(meanVulPrey,'omitnan');
meanVulPreyFree = wmean(meanVulPrey,free)/mean(meanVulPrey(binCons));
meanVulPreyPara = wmean(meanVulPrey,para)/mean(meanVulPrey(binCons));


%PreyAveraged Trophic LEvel
A(eye(S)>0) = 0;
A = sparse(A);
paTL_mx = sparse(A*(diag(1./sum(A))));
patl = (speye(S)-paTL_mx')\ones(S,1);

patlFree = wmean(patl,free)/mean(patl(binCons));
patlPara = wmean(patl,para)/mean(patl(binCons));

%What follows is various ccs from an arxiv paper:
%Clustering Coefficient(s)
A = A*1;
A2 = full(A^2);
d2 = diag(A2);
nNayb = gen0+vul0;
Asym = A+A';
cc0 = full(diag((Asym)^3)./(2*(nNayb).*(nNayb-1)-2*d2));
cc0(isnan(cc0))= 0;
cc0Free = wmean(cc0,free)/mean(cc0(binCons));
cc0Para = wmean(cc0,para)/mean(cc0(binCons));

%Cyclic; 
ccCyc = full(diag(A*A2))./(gen0.*vul0-d2);
ccCyc(isnan(ccCyc)) = 0;
ccCycFree = wmean(ccCyc,free)/mean(ccCyc(binCons));
ccCycPara = wmean(ccCyc,para)/mean(ccCyc(binCons));
%middleman
ccMid = full(diag((A*A')*A))./(gen0.*vul0-d2);
ccMid(isnan(ccMid)) = 0;
ccMidFree = wmean(ccMid,free)/mean(ccMid(binCons));
ccMidPara = wmean(ccMid,para)/mean(ccMid(binCons));
%innie
ccIn = full(diag(A'*A2))./(gen0.*(gen0-1));
ccIn(isnan(ccIn)) = 0;
ccInFree = wmean(ccIn,free)/mean(ccIn(binCons));
ccInPara = wmean(ccIn,para)/mean(ccIn(binCons));
%outie
ccOut = full(diag(A2*A'))./(vul0.*(vul0-1));
ccOut(isnan(ccOut)) = 0;
ccOutFree = wmean(ccOut,free)/mean(ccOut(binCons));
ccOutPara = wmean(ccOut,para)/mean(ccOut(binCons));

%betweenness takes a ton of time; would be better if I could do this 
%faster.. oh well.

%Betweenness
btwns = btwn(res,con);
btwns = (btwns-min(btwns))/max(btwns);
btwnsFree = wmean(btwns,free)/mean(btwns(binCons));
btwnsPara = wmean(btwns,para)/mean(btwns(binCons));

%Ecological Betweenness
ecoBtwns = ecoBtwn(res,con);
ecoBtwns = (ecoBtwns-min(ecoBtwns))/max(ecoBtwns);
ecoBtwnsFree = wmean(ecoBtwns,free)/mean(ecoBtwns(binCons));
ecoBtwnsPara = wmean(ecoBtwns,para)/mean(ecoBtwns(binCons));

%Ecological Pagerank
pr = pageRank(res,con);
prFree = wmean(pr,free)/mean(pr(binCons));
prPara = wmean(pr,para)/mean(pr(binCons));

C = L/S^2;

Lcon = sum(sum(A(binCons,binCons)));
Scon = sum(binCons);
Ccon = Lcon/Scon^2;

binPara = para>=0.5;
binFree = free>0.5;

Lff = sum(sum(A.*(free*free')));
Sf = sum(binFree);
Cff = Lff/Sf^2;

Lpp = sum(sum(A.*(para*para')));
Sp = sum(binPara);
Cpp = Lpp/Sp^2;

Lfp = sum(sum(A.*(free*para')));
Cfp = Lfp/(Sp*Sf);

Lpf = sum(sum(A.*(para*free')));
Cpf = Lpf/(Sp*Sf);

fPar = Sp/(Sp+Sf);

top = sum(vul == 0)/S;
bas = sum(gen==0)/S;
int = 1-top-bas;
herb = sum(abs(patl-2)<1e-2)/S;
omn = sum(mod(patl,1)>1e-2)/S;

para = para./nodeSize;
free = free./nodeSize;

localProperties = [     
                         vul...                 1
                        ,gen...                 2
                        ,patl...                3
                        ,cc0...                 4
                        ,pr...                  5
                        ,btwns...               6
                        ,ecoBtwns...            7
                        ,meanVulPrey...         8
                        ,meanGenPred...         9
                        ,para...                10
                        ,free...                11
                        ,repmat(S,S,1)...       12
                        ,repmat(C,S,1)...       13
                        ,nodeSize...            14
                        ,ccCyc...                 15
                        ,ccMid...                 16
                        ,ccIn...                 17
                        ,ccOut...                 18
                    ];

globalProperties = [
                   S...             1   
                  ,Ccon...          2
                  ,Cff...           3
                  ,Cpp...           4
                  ,Cfp...           5
                  ,Cpf...           6
                  ,fPar...          7
                  ,numLargest...    8  
                  ,n...             9
                  ,top...           10
                  ,int...           11
                  ,bas...           12
                  ,herb...          13
                  ,omn...           14
                  ,C...             15
                  ,minDistance...   16
                      ];     
localMeans = [                 
                 vulFree...                 1
                ,vulPara...                 2
                ,genFree...                 3
                ,genPara...                 4
                ,patlFree...                5
                ,patlPara...                6
                ,cc0Free...                  7
                ,cc0Para...                  8
                ,prFree...                  9
                ,prPara...                  10
                ,btwnsFree...               11
                ,btwnsPara...               12
                ,ecoBtwnsFree...            13
                ,ecoBtwnsPara...            14
                ,meanVulPreyFree...         15
                ,meanVulPreyPara...         16
                ,meanGenPredFree...         17
                ,meanGenPredPara...         18
                ,S...                       19
                ,C...                       20
                ,minDistance...             21
                ,ccCycFree...                 22
                ,ccCycPara...                 23
                ,ccMidFree...                 24
                ,ccMidPara...                 25
                ,ccInFree...                 26
                ,ccInPara...                 27
                ,ccOutFree...                 28
                ,ccOutPara...                 29
            ];
end