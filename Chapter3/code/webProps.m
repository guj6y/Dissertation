function [localProperties,globalProperties,localMeans] = webProps(A,para)

%Now, we calculate properties using the largest connected component.
S = length(A);
L = sum(sum(A));

%In degree (generality)
gen0 = full(sum(A))';

basal = gen0==0;
binFree = ~binPara;
binFree(basal) = false;
binCons = (~basal);
%Identify carnivores Only!
carn = true(S,1);
carn(con(basal(res))) = false;
carn(para) = false;

gen = (gen0)./mean(gen0);
genFree = mean(gen(carn)));
genPara = mean(gen(para)));

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
meanGenPredFree = mean(meanGenPred(carn));
meanGenPredPara = mean(meanGenPred(para));

%out degree (vulnerability)
vul0 = full(sum(A,2));
vul = (vul0)/mean(vul0);
vulFree = mean(vul(carn));
vulPara = mean(vul(para));


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
meanVulPreyFree = mean(meanVulPrey(carn));
meanVulPreyPara = mean(meanVulPrey(para));


%PreyAveraged Trophic LEvel
A(eye(S)>0) = 0;
A = sparse(A);
paTL_mx = sparse(A*(diag(1./sum(A))));
patl = (speye(S)-paTL_mx')\ones(S,1);
patl = (patl-1)/mean(patl-1);
patlFree = mean(patl(carn));
patlPara = mean(patl(para));



%betweenness takes a ton of time; would be better if I could do this 
%faster.. oh well.

%Betweenness
btwns = btwn(res,con);
btwns = (btwns)/mean(btwns,'omitnan');
btwnsFree = mean(btwns(carn));
btwnsPara = mean(btwns(para));

%Ecological Betweenness
ecoBtwns = ecoBtwn(res,con);
ecoBtwns = (ecoBtwns)/mean(ecoBtwns,'omitnan');
ecoBtwnsFree = mean(ecoBtwns(carn));
ecoBtwnsPara = mean(ecoBtwns(para));

%Ecological Pagerank
pr = log(pageRank(res,con));
pr = (pr)/mean(pr,'omitnan');
prFree = mean(pr(carn));
prPara = mean(pr(para));

%What follows is various ccs from an arxiv paper:
%Clustering Coefficient(s)
A = A*1;
A2 = full(A^2);
d2 = diag(A2);
nNayb = gen0+vul0;
Asym = A+A';
cc0 = full(diag((Asym)^3)./(2*(nNayb).*(nNayb-1)-2*d2));

%cc0(isnan(cc0))= 0;
cc0Free = mean(cc0(carn))/mean(cc0(binCons),'omitnan');
cc0Para = mean(cc0(para))/mean(cc0(binCons),'omitnan');

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
para = para + 0;
binCons = binCons+0;

Lff = sum(sum(A.*(binCons*binCons')));
Sf = sum(binFree);
Cff = Lff/Sf^2;

Lpp = sum(sum(A.*(para*para')));
Sp = sum(binPara);
Cpp = Lpp/Sp^2;

Lfp = sum(sum(A.*(binCons*para')));
Cfp = Lfp/(Sp*Sf);

Lpf = sum(sum(A.*(para*binCons')));
Cpf = Lpf/(Sp*Sf);

fPar = Sp/(Sp+Sf);

top = sum(vul == 0)/S;
bas = sum(gen==0)/S;
int = 1-top-bas;
herb = sum(abs(patl-2)<1e-2)/S;
omn = sum(mod(patl,1)>1e-2)/S;

para = para./nodeSize;
carn = carn./nodeSize;

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
                      ];     

end