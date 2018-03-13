function [localProperties,globalProperties,pearCorr,spearCorr] = agglomProps(A,para)

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
    localProperties = zeros(1,18);
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
gen = full(sum(A));
gen = gen';

%This defines parasite cutoffs; without this line we calculate weighted
%averages:
binPara = para>0.5;
basal = gen==0;
binFree = ~binPara;
binFree(basal) = false;
binCons = (~basal);

free = 1-para;
free(basal) = 0;

gen = gen./(mean(gen(binCons)));

genPara = wmean(gen,para);
genFree = wmean(gen,free);

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
meanGenPred = meanGenPred./mean(meanGenPred(binCons));
meanGenPredFree = wmean(meanGenPred,free);
meanGenPredPara = wmean(meanGenPred,para);

%out degree (vulnerability)
vul = full(sum(A,2));
vul = vul./mean(vul(binCons));

vulPara = wmean(vul,para);
vulFree = wmean(vul,free);

%Mean vulnerability of prey
%Mean topological fraction of preys' binConsumers
meanVulPrey = zeros(S,1);
for ii = 1:S
    prey_ii = res(con==ii);
    vulPrey_ii = vul(prey_ii);
    meanVulPrey(ii) = mean(vulPrey_ii);
end

meanVulPrey = meanVulPrey/mean(meanVulPrey(binCons));
meanVulPreyFree = wmean(meanVulPrey,free);
meanVulPreyPara = wmean(meanVulPrey,para);

if sum(isnan(meanVulPrey)>0)
    fprintf('wtfmate\n')
    
end

%Formula for the prey-averaged trophic level.
%no cannibalism in this matrix; think it's better this way:
A(((1:S)-1)*S+(1:S)) = 0; 
patl_mx = sparse(A)*(diag(1./sum(sparse(A))));

patl = (speye(S)-patl_mx')\ones(S,1);
patl = patl./mean(patl(binCons));

patlPara = wmean(patl,para);
patlFree = wmean(patl,free);

%Clustering Coefficient(s)
cc = clustering_coefficients(sparse(A));
cc = cc./mean(cc(binCons));

ccPara = wmean(cc,para);
ccFree = wmean(cc,free);

%Betweenness
btwns = btwn(res,con);
btwns = btwns./mean(btwns(binCons));
btwnFree = wmean(btwns,free);
btwnPara = wmean(btwns,para);

%Ecological Betweenness
ecoBtwns = ecoBtwn(res,con);
ecoBtwns = ecoBtwns./mean(ecoBtwns(binCons));
ecoBtwnFree = wmean(ecoBtwns,free);
ecoBtwnPara = wmean(ecoBtwns,para);

%Ecological Pagerank
pr = pageRank(res,con);
pr = pr./mean(pr(binCons));
prFree = wmean(pr,free);
prPara = wmean(pr,para);


localProperties = [     
                         vulFree...                 1
                        ,vulPara...                 2
                        ,genFree...                 3
                        ,genPara...                 4
                        ,patlFree...                5
                        ,patlPara...                6
                        ,ccFree...                  7
                        ,ccPara...                  8
                        ,prFree...                  9
                        ,prPara...                  10
                        ,btwnFree...                11
                        ,btwnPara...                12
                        ,ecoBtwnFree...             13
                        ,ecoBtwnPara...             14
                        ,meanVulPreyFree...         15
                        ,meanVulPreyPara...         16
                        ,meanGenPredFree...         17
                        ,meanGenPredPara...         18
                    ];

pearCorr = corr([vul,gen,patl,cc,pr,btwns,ecoBtwns,meanVulPrey,meanGenPred]);
spearCorr = corr([vul,gen,patl,cc,pr,btwns,ecoBtwns,meanVulPrey,meanGenPred],'type','Spearman');
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
herb = sum(patl==2)/S;
omn = sum(mod(patl,1)>.01)/S;

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