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

%eco btwons and pagerank are the rate limiting step.
ecoBtwns = ecoBtwn(res,con);

%Ecological Pagerank
pr = (pageRank(res,con));

X=  [vul...                 1
    ,gen...                 2
    ,patl...                3
    ,pr...                  4
    ,meanVulPrey...         5
    ,meanGenPred...         6
    ,ecoBtwns...            7
    ];

X(isnan(X)) = 0;

end
