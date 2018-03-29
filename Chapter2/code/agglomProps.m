function [localProperties,globalProperties,localMeans,res,con] = agglomProps(A,para,nodeSize,minDistance,n)

%Check for the connected component; don't worry about the web having basal
%species, etc. since we are not trying to develop a dynamical model for
%these; we would deal with that in the modeling phase anyway. We do want to
%get rid of disconnected elements, though. This should (hopefully) have the
%added benefit that we don't get as many warnings when calculating the
%trophic leve.
S0 = length(A);

[n,grps] = graphconncomp(sparse(A),'weak',true,'directed',true);

compSizes = zeros(S0,1);
nLinks = zeros(S0,1);
for ii = 1:n
    thisGrp = grps==ii;
    sizeii = sum(thisGrp);
    nLinks(ii) = sum(sum(A(thisGrp,thisGrp)));
    compSizes(ii) = sizeii;
end

sortedComps = sortrows([(1:S0)', compSizes, nLinks, rand(S0,1)],[-2 -3 4]);
numLargest = sum(sortedComps(:,2) == sortedComps(1,2));

winners = grps==sortedComps(1);
%I like winners. Winners get to stay.

A = A(winners,winners);
[res,con] = find(A);

if isempty(res) == 1
    localProperties = zeros(1,19);
    globalProperties = zeros(1,20);
    localMeans = zeros(1,30);
    res = 0;
    con = 0;
    return
end
para = para(winners);
%carn shouldn't be carried around like this; it is a STRUCTURAL property!
%(though it's weird maybe that parasite aren't staying carnivores..?)
%carn = carn(winners);
nodeSize = nodeSize(winners);



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
binCons = (~basal);
%Identify carnivores Only!
carn = true(S,1);
carn(con(basal(res))) = false;
carn(basal) = false;
carn(binPara) = false;


weighted=false;
if weighted
    %Weigthing Option:
    %Pick the appropriate other comparison.
    %free = free*nodeSize;
    carn = carn.*nodeSize;
    para = para.*nodeSize;
else
    %Majority Option:
    carn = carn>0.5;
    para = para>=0.5;
end
gen = (gen0)./mean(gen0);
genFree = wmean(gen,carn);
genPara = wmean(gen,para);

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
meanGenPred = (meanGenPred)/mean(meanGenPred,'omitnan');
meanGenPredFree = wmean(meanGenPred,carn);
meanGenPredPara = wmean(meanGenPred,para);

%out degree (vulnerability)
vul0 = full(sum(A,2));
vul = (vul0)/mean(vul0);
vulFree = wmean(vul,carn);
vulPara = wmean(vul,para);


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

meanVulPrey = (meanVulPrey-mean(meanVulPrey,'omitnan'))/std(meanVulPrey,'omitnan');
meanVulPreyFree = wmean(meanVulPrey,carn);
meanVulPreyPara = wmean(meanVulPrey,para);


%PreyAveraged Trophic LEvel
A(eye(S)>0) = 0;
A = sparse(A);
paTL_mx = sparse(A*(diag(1./sum(A))));
patl = (speye(S)-paTL_mx')\ones(S,1);
patl = (patl-1)/mean(patl-1);
patlFree = wmean(patl,carn);
patlPara = wmean(patl,para);



%betweenness takes a ton of time; would be better if I could do this 
%faster.. oh well.

%Betweenness
btwns = btwn(res,con);
btwns = (btwns)/mean(btwns,'omitnan');
btwnsFree = wmean(btwns,carn);
btwnsPara = wmean(btwns,para);

%Ecological Betweenness
ecoBtwns = ecoBtwn(res,con);
ecoBtwns = (ecoBtwns)/mean(ecoBtwns,'omitnan');
ecoBtwnsFree = wmean(ecoBtwns,carn);
ecoBtwnsPara = wmean(ecoBtwns,para);

%Ecological Pagerank
pr = log(pageRank(res,con));
pr = (pr)/mean(pr,'omitnan');
prFree = wmean(pr,carn);
prPara = wmean(pr,para);

%What follows is various ccs from an arxiv paper:
%Clustering Coefficient(s)
A = A*1;
A2 = full(A^2);
d2 = diag(A2);
nNayb = gen0+vul0;
Asym = A+A';
cc0 = full(diag((Asym)^3)./(2*(nNayb).*(nNayb-1)-2*d2));

%cc0(isnan(cc0))= 0;
cc0Free = wmean(cc0,carn)/mean(cc0(binCons),'omitnan');
cc0Para = wmean(cc0,para)/mean(cc0(binCons),'omitnan');

%Cyclic; 
ccCycTop = full(diag(A*A2));
ccCycBot = (gen0.*vul0-d2);
ccCyc = ccCycTop./ccCycBot;
ccCyc = ccCyc/mean(ccCyc,'omitnan');
%ccCyc(isnan(ccCyc)) = 0;
ccCycCarn = sum(ccCycTop.*carn)/sum(ccCycBot.*carn);
ccCycPara = sum(ccCycTop.*para)/sum(ccCycBot.*para);


%middleman
ccMidTop = full(diag((A*A')*A));
ccMidBot = (gen0.*vul0-d2);
ccMid = ccMidTop./ccMidBot;
ccMid = ccMid/mean(ccMid,'omitnan');
%ccMid(isnan(ccMid)) = 0;
ccMidCarn = sum(ccMidTop.*carn)/sum(ccMidBot.*carn);
ccMidPara = sum(ccMidTop.*para)/sum(ccMidBot.*para);


%innie
ccInTop = full(diag(A'*A2));
ccInBot = (gen0.*(gen0-1));
ccIn = ccInTop./ccInBot;
ccIn = ccIn/mean(ccIn,'omitnan');

%ccIn(isnan(ccIn)) = 0;
ccInCarn = sum(ccInTop.*carn)/sum(ccInBot.*carn);
ccInPara = sum(ccInTop.*para)/sum(ccInBot.*para);
%outie
ccOutTop = full(diag(A2*A'));
ccOutBot = (vul0.*(vul0-1));
ccOut = ccOutTop./ccOutBot;
ccOut = ccOut/mean(ccOut,'omitnan');

%ccOut(isnan(ccOut)) = 0;
ccOutCarn = sum(ccOutTop.*carn)/sum(ccOutBot.*carn);
ccOutPara = sum(ccOutTop.*para)/sum(ccOutBot.*para);

C = L/S^2;

Lcon = sum(sum(A(binCons,binCons),'omitnan'));
Scon = sum(binCons);
Ccon = Lcon/Scon^2;

carn = carn + 0;
free = ~para;

Lff = sum(free(res)&free(con));
Sf = sum(free);
Cff = Lff/Sf^2;

Lpp = sum(para(res)&para(con));
Sp = sum(para);
Cpp = Lpp/Sp^2;

Lfp = sum(free(res)&para(con));
Cfp = Lfp/(Sp*Sf);

Lpf = sum(para(res)&free(con));
Cpf = Lpf/(Sp*Sf);

Cf = (Lff + Lpf)/(Sf*(Sp+Sf));
Cp = (Lpp + Lfp)/(Sp*(Sp+Sf));

fPar = Sp/(Sp+Sf);


top = sum(vul == 0)/S;
bas = mean(basal);
int = 1-top-bas;
haveBasalRes = false(S,1);
haveBasalRes(con(basal(res))) = true;
haveConRes = false(S,1);
haveConRes(con(~basal(res))) = true;
herbs = haveBasalRes&(~haveConRes);
carns = haveConRes&(~haveBasalRes);
omns = haveBasalRes&haveConRes;
herb = mean(herbs);
omn = mean(omns);

can = sum(res==con)/S;
fCarn = mean(carn);


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
                        ,carn...                11
                        ,repmat(S,S,1)...       12
                        ,repmat(C,S,1)...       13
                        ,nodeSize...            14
                        ,ccCyc...               15
                        ,ccMid...               16
                        ,ccIn...                17
                        ,ccOut...               18
                        ,repmat(n,S,1)...       19
                    ];
                


globalProperties = [
                   S...             1   
                  ,Cf...            2
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
                  ,fCarn...         17
                  ,Cp...            18
                  ,Sf...            19
                  ,n...             20
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
                ,ccCycCarn...                 22
                ,ccCycPara...                 23
                ,ccMidCarn...                 24
                ,ccMidPara...                 25
                ,ccInCarn...                 26
                ,ccInPara...                 27
                ,ccOutCarn...                 28
                ,ccOutPara...                 29
                ,n...                         30
            ];
end