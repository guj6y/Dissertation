function [X] = propsForClassTree(res,con)


S = max([res;con]);
%Now, we calculate properties using the largest connected component.
A = sparse(res,con,1,S,S);
%In degree (generality)
gen = full(sum(A))';


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


%out degree (vulnerability)
vul = full(sum(A,2));



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



%PreyAveraged Trophic LEvel
A(eye(S)>0) = 0;
A = sparse(A);
paTL_mx = sparse(A*(diag(1./sum(A))));
patl = (speye(S)-paTL_mx')\ones(S,1);


%betweenness takes a ton of time; would be better if I could do this 
%faster.. oh well.

%Betweenness
btwns = betweenness_centrality(A);

%Ecological Betweenness
ecoBtwns = ecoBtwn(res,con);

%Ecological Pagerank
pr = log(pageRank(res,con));


%What follows is various ccs from an arxiv paper:
%Clustering Coefficient(s)
A = A*1;
A2 = full(A^2);
d2 = diag(A2);
%Cyclic; 
ccCycTop = full(diag(A*A2));
ccCycBot = (gen.*vul-d2);
ccCyc = ccCycTop./ccCycBot;

%middleman
ccMidTop = full(diag((A*A')*A));
ccMidBot = (gen.*vul-d2);
ccMid = ccMidTop./ccMidBot;

%innie
ccInTop = full(diag(A'*A2));
ccInBot = (gen.*(gen-1));
ccIn = ccInTop./ccInBot;

%outie
ccOutTop = full(diag(A2*A'));
ccOutBot = (vul.*(vul-1));
ccOut = ccOutTop./ccOutBot;

X=  [vul...                 1
    ,gen...                 2
    ,patl...                3
    ,pr...                  4
    ,btwns...               5
    ,ecoBtwns...            6
    ,meanVulPrey...         7
    ,meanGenPred...         8
    ,ccCyc...               9
    ,ccMid...               10
    ,ccIn...                11
    ,ccOut...               12
    ];

end