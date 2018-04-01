function [output] = webPropsRaw(res,con,para)
%non-normalized properties!
%Now, we calculate properties using the largest connected component.
S = max([res;con]);

A = sparse(res,con,1,S,S);
X = propsForClassTree(res,con);

X(:,[1:2 4:12]) = X(:,[1:2 4:12])./mean(X(:,[1:2 4:12]),'omitnan');
X(:,3) = (X(:,3)-1)/mean(X(:,3)-1);
X(isnan(X)) = 0;
vul = X(:,1);
gen = X(:,2);
patl = X(:,3);
pr = X(:,4);
btwns = X(:,5); 
ecoBtwns = X(:,6);
meanVulPrey = X(:,7);
meanGenPred = X(:,8);

basal = gen==0;
%Identify carnivores Only!
carn = true(S,1);
carn(con(basal(res))) = false;
carn(basal) = false;
carn(para) = false;

genFree = mean(gen(carn));
genPara = mean(gen(para));
genDiff = genFree-genPara;

meanGenPredFree = mean(meanGenPred(carn));
meanGenPredPara = mean(meanGenPred(para));
meanGenPredDiff = meanGenPredFree-meanGenPredPara;

%out degree (vulnerability)
vulFree = mean(vul(carn));
vulPara = mean(vul(para));
vulDiff = vulFree-vulPara;


%Mean vulnerability of prey
%Mean topological fraction of preys' binConsumers
meanVulPreyFree = mean(meanVulPrey(carn));
meanVulPreyPara = mean(meanVulPrey(para));
meanVulPreyDiff = meanVulPreyFree - meanVulPreyPara;

%PreyAveraged Trophic LEvel
patlFree = mean(patl(carn));
patlPara = mean(patl(para));
patlDiff = patlFree - patlPara;


%betweenness takes a ton of time; would be better if I could do this 
%faster.. oh well.

%Betweenness
btwnsFree = mean(btwns(carn));
btwnsPara = mean(btwns(para));
btwnsDiff = btwnsFree-btwnsPara;

%Ecological Betweenness
ecoBtwnsFree = mean(ecoBtwns(carn));
ecoBtwnsPara = mean(ecoBtwns(para));
ecoBtwnsDiff = ecoBtwnsFree-ecoBtwnsPara;

%Ecological Pagerank
prFree = mean(pr(carn));
prPara = mean(pr(para));
prDiff = prFree-prPara;

A = A*1;
A(eye(S)>0) = 0;
A2 = full(A^2);
d2 = diag(A2);

%Cyclic; 
ccCycTop = full(diag(A*A2));
ccCycBot = (gen.*vul-d2);
ccCyc = sum(ccCycTop)./sum(ccCycBot);
%ccCyc(isnan(ccCyc)) = 0;
ccCycCarn = sum(ccCycTop.*carn)/sum(ccCycBot.*carn);
ccCycPara = sum(ccCycTop.*para)/sum(ccCycBot.*para);
ccCycDiff = ccCycCarn-ccCycPara;

%middleman
ccMidTop = full(diag((A*A')*A));
ccMidBot = (gen.*vul-d2);
ccMid = sum(ccMidTop)./sum(ccMidBot);
%ccMid(isnan(ccMid)) = 0;
ccMidCarn = sum(ccMidTop.*carn)/sum(ccMidBot.*carn);
ccMidPara = sum(ccMidTop.*para)/sum(ccMidBot.*para);
ccMidDiff = ccMidCarn-ccMidPara;

%innie
ccInTop = full(diag(A'*A2));
ccInBot = (gen.*(gen-1));
ccIn = sum(ccInTop)./sum(ccInBot);

%ccIn(isnan(ccIn)) = 0;
ccInCarn = sum(ccInTop.*carn)/sum(ccInBot.*carn);
ccInPara = sum(ccInTop.*para)/sum(ccInBot.*para);
ccInDiff = ccInCarn-ccInPara;

%outie
ccOutTop = full(diag(A2*A'));
ccOutBot = (vul.*(vul-1));
ccOut = sum(ccOutTop)./sum(ccOutBot);

%ccOut(isnan(ccOut)) = 0;
ccOutCarn = sum(ccOutTop.*carn)/sum(ccOutBot.*carn);
ccOutPara = sum(ccOutTop.*para)/sum(ccOutBot.*para);
ccOutDiff = ccOutCarn-ccOutPara;

top = sum(vul == 0)/S;
bas = mean(basal);
int = 1-top-bas;
herbs = false(S,1);
herbs(con(basal(res))) = true;
herbs(con(~basal(res))) = false;
herb = mean(herbs);
omn = mean(~(basal|herbs|carn));
can = sum(res==con)/S;
fCarn = mean(carn);

linkSd = std(gen+vul);
genSd = std(gen);
vulSd = std(vul);
TL = mean(patl);
loop = 0;
[~,grps] = graphconncomp(sparse(A),'directed',true);
for ii = 1:S
    if sum(grps==grps(ii))>1
        loop = loop+1;
    end
    
end
loop = loop/S;
maxSim = mean(max(calculateSimilarity(res,con)));
A_ = ((A+A'));
A_(A_>0) = 1.0;

D = all_shortest_paths(A_);
path = mean(D(:));

globalProps =      [top...      1
                   ,int...      2
                   ,bas...      3
                   ,herb...     4
                   ,omn...      5
                   ,can...      6
                   ,loop...     7
                   ,linkSd...   8
                   ,genSd...    9
                   ,vulSd...    10
                   ,TL...       11
                   ,maxSim...   12
                   ,path...     13
                   ,ccCyc...    14
                   ,ccMid...    15
                   ,ccIn...     16
                   ,ccOut...    17
                   ,fCarn...     18
                   ];
                      
localDiffs = [   vulDiff...                 1
                ,genDiff...                 2
                ,patlDiff...                3
                ,prDiff...                  4
                ,btwnsDiff...               5
                ,ecoBtwnsDiff...            6
                ,meanVulPreyDiff...         7
                ,meanGenPredDiff...         8
                ,ccCycDiff...               9
                ,ccMidDiff...               10
                ,ccInDiff...                11
                ,ccOutDiff...               12
            ];
globalProps(isnan(globalProps)) = 0;
localDiffs(isnan(localDiffs)) = 0;

output = {globalProps, localDiffs, X};
end