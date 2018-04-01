function [output] = webPropsRaw(res,con,para)

%non-normalized properties!
%Now, we calculate properties using the largest connected component.
S = max([res;con]);

A = sparse(res,con,1,S,S);
X = propsForClassTree(res,con);

vul = X(:,1);
gen = X(:,2);
patl = X(:,3);
pr = X(:,4);
ecoBtwns = X(:,7);
meanVulPrey = X(:,5);
meanGenPred = X(:,6);

basal = gen==0;

free=~para;

genFree = mean(gen(free));
genPara = mean(gen(para));

meanGenPredFree = mean(meanGenPred(free));
meanGenPredPara = mean(meanGenPred(para));

%out degree (vulnerability)
vulFree = mean(vul(free));
vulPara = mean(vul(para));

%Mean vulnerability of prey
%Mean topological fraction of preys' binConsumers
meanVulPreyFree = mean(meanVulPrey(free));
meanVulPreyPara = mean(meanVulPrey(para));

%betweenness takes a ton of time; would be better if I could do this 
%faster.. oh well.

%Ecological Betweenness
ecoBtwnsFree = mean(ecoBtwns(free));
ecoBtwnsPara = mean(ecoBtwns(para));
ecoBtwnsAll = mean(ecoBtwns);

%Ecological Pagerank
prFree = mean(log10(pr(free)));
prPara = mean(log10(pr(para)));

prSD = std(log10(pr));

A = A*1;
A(eye(S)>0) = 0;
A2 = full(A^2);
d2 = diag(A2);

%Cyclic; 
ccCycTop = full(diag(A*A2));
ccCycBot = (gen.*vul-d2);
ccCyc = sum(ccCycTop)./sum(ccCycBot);
%ccCyc(isnan(ccCyc)) = 0;
ccCycFree = sum(ccCycTop.*free)/sum(ccCycBot.*free);
ccCycPara = sum(ccCycTop.*para)/sum(ccCycBot.*para);

%middleman
ccMidTop = full(diag((A*A')*A));
ccMidBot = (gen.*vul-d2);
ccMid = sum(ccMidTop)./sum(ccMidBot);
%ccMid(isnan(ccMid)) = 0;
ccMidFree = sum(ccMidTop.*free)/sum(ccMidBot.*free);
ccMidPara = sum(ccMidTop.*para)/sum(ccMidBot.*para);

%innie
ccInTop = full(diag(A'*A2));
ccInBot = (gen.*(gen-1));
ccIn = sum(ccInTop)./sum(ccInBot);
%ccIn(isnan(ccIn)) = 0;
ccInFree = sum(ccInTop.*free)/sum(ccInBot.*free);
ccInPara = sum(ccInTop.*para)/sum(ccInBot.*para);

%outie
ccOutTop = full(diag(A2*A'));
ccOutBot = (vul.*(vul-1));
ccOut = sum(ccOutTop)./sum(ccOutBot);
%ccOut(isnan(ccOut)) = 0;
ccOutFree = sum(ccOutTop.*free)/sum(ccOutBot.*free);
ccOutPara = sum(ccOutTop.*para)/sum(ccOutBot.*para);

%Identify carnivores Only!
carn = true(S,1);
carn(con(basal(res))) = false;
carn(basal) = false;
carn(para) = false;

tops = vul == 0;
top = mean(tops);
topsFree = mean(tops(free));
topsPara = mean(tops(para));

bas = mean(basal);
int = 1-top-bas;
intFree = 1-topFree - mean(basal(free));
intPara = 1-topPara; 

herbs = false(S,1);
herbs(con(basal(res))) = true;
herbs(con(~basal(res))) = false;
herb = mean(herbs);
herbFree = mean(herbs(free));
herbPara = mean(herbs(para));

omns = ~(basal|herbs|carn);
omn = mean(omns);
omnFree = mean(omns(free));
omnPara = mean(omns(para));

can = sum(res==con)/S;

fCarn = mean(carn);
fCarnFree = mean(carn(free));
fCarnPara = mean(carn(para));

genSd = std(gen);
vulSd = std(vul);
rgv = corr(gen,vul);

genSdFree = std(gen(free));
vulSdFree = std(vul(free));
rgvFree = corr(gen(free),vul(free));

genSdPara = std(gen(para));
vulSdPara = std(vul(para))
rgvFree = corr(gen(para),vul(para));


TL = mean(patl);
TLFree = mean(patl(free));
TLPara = mean(patl(para));

loop = 0;
loopFree = 0;
loopPara = 0;
[~,grps] = graphconncomp(sparse(A),'directed',true);

for ii = 1:S
    if sum(grps==grps(ii))>1
        loop = loop+1;
        loopFree = loop + free(ii);
        loopPara = loop + para(ii);
    end
    
end
loop = loop/S;
loopFree = loopFree/sum(free);
loopPara = loopPara/sum(para);

simMx = calculateSimilarity(res,con))
maxSim = mean(max(simMx);
maxSimFree = mean(max(simMx(free,free)));
maxSimPara = mean(max(simMx(para,para)));

A_ = ((A+A'));
A_(A_>0) = 1.0;
D = all_shortest_paths(A_);
D = triu(D,1);
DFree = D(free,free);
DPara = D(para,para);

path = mean(D(D>0));
pathFree = mean(DFree(DFree>0));
pathPara = mean(DPara(DPara>0));



globalProps =      [top...          1
                   ,int...          2
                   ,bas...          3
                   ,herb...         4
                   ,omn...          5
                   ,fCarn...        6
                   ,can...          7
                   ,rgv...          8
                   ,genSd...        9
                   ,vulSd...        10
                   ,TL...           11
                   ,prSD...         12 
                   ,path...         13
                   ,loop...         14
                   ,ecoBtwnsAll...  15
                   ,ccCyc...        16
                   ,ccMid...        17
                   ,ccIn...         18
                   ,ccOut...        19
                   ,maxSim...       20
                   ];
                      
freeProps =  [   vulFree...                 1
                ,genFree...                 2
                ,TLFree...                  3
                ,prFree...                  4
                ,ecoBtwnsFree...            5
                ,meanVulPreyFree...         6
                ,meanGenPredFree...         7
                ,ccCycFree...               8
                ,ccMidFree...               9
                ,ccInFree...                10
                ,ccOutFree...               11
                ,topsFree...                12
                ,intFree...                 13
                ,herbFree...                14
                ,omnFree...                 15
                ,fCarnFree...               16
                ,genSdFree...               17
                ,vulSdFree...               18
                ,rgvFree...                 19
                ,loopFree...                20
                ,maxSimPara...              21
                ,pathFree...                22

            ];

paraProps = [    vulPara...                 1
                ,genPara...                 2
                ,TLPara...                  3
                ,prPara...                  4
                ,ecoBtwnsPara...            5
                ,meanVulPreyPara...         6
                ,meanGenPredPara...         7
                ,ccCycPara...               8
                ,ccMidPara...               9
                ,ccInPara...                10
                ,ccOutPara...               11
                ,topsPara...                12
                ,intPara...                 13
                ,herbPara...                14
                ,omnPara...                 15
                ,fCarnPara...               16
                ,genSdPara...               17
                ,vulSDPara...               18
                ,rgvPara...                 19
                ,loopPara...                20
                ,maxSimPara...              21
                ,pathPara...                22
            ];


globalProps(isnan(globalProps)) = 0;
freeProps(isnan(freeProps)) = 0;
paraProps(isnan(paraProps)) = 0;

output = {globalProps, freeProps, paraProps};

end
