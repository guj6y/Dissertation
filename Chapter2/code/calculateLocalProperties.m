function [local_properties] = calculateLocalProperties(res,con)

N = max(max(res,con));
A = sparse(res,con,1.0,N,N);
basalsp = full(sum(A)==0);
nonbasalsp = ~basalsp;
nicheweb = full(A);
S = 1:N;
L = sum(sum(nicheweb));
%Clustering Coefficient(s)
cc = clustering_coefficients(A);

%In degree (generalit)
gen = full(sum(A))/(L/N);

%out degree (vulnerability)
vul = full(sum(A,2))/(L/N);

%Mean vulnerability of prey
%Mean topological fraction of preys' consumers
mean_vul_prey = zeros(N,1);

for ii = 1:N
    
    prey_ii = nicheweb(:,ii)>0;
    vul_prey_ii = sum(nicheweb(prey_ii,:),2)/(L/N);
    mean_vul_prey(ii) = mean(vul_prey_ii);
    
end

%Mean generality of predators
%Mean topological fraction of predators' diet
mean_gen_pred = zeros(N,1);
for ii = 1:N
    
    pred_ii = nicheweb(ii,:)>0;
    gen_pred_ii = sum(nicheweb(:,pred_ii))/(L/N);
    mean_gen_pred(ii) = mean(gen_pred_ii);
end

[betweenness,~,~,min_sp_tobasal,numBasalCon] = calcEcologicalBetweenness(res,con);

frac_basal_con = numBasalCon./sum(basalsp);

%short-weighted trophic level

%Formula for the prey-averaged trophic level.
paTL_mx = A*(diag(1./sum(A)));

paTL = (speye(N)-paTL_mx')\ones(N,1);

%Short trophic level is shortest path to a basal; since I've already
%calculated shortest paths this is easy >;* <<<3 xoxo -nk
sTL = min_sp_tobasal+1;

swTL = (paTL+sTL)/2;

%Number of species in loops is the number of species in strongly connected
%components larger than one species.
[n_comp,comps] = graphconncomp(A,'Directed',1,'Weak',0);
num_comp = zeros(1,n_comp);
sp_in_loop = [];

for ii = 1:n_comp
    num_comp(ii) = sum(comps==ii);
    if num_comp(ii) > 1
        sp_in_loop = [sp_in_loop; S(comps==ii)'];
    end
end
inLoop = zeros(N,1);
inLoop(sp_in_loop) = 1;

top_sp = sum(A,2)==0;
%basal?
mean_vul_prey(basalsp) = 0;

mean_gen_pred(top_sp) = 0;
pageRank = calculatePageRank(res,con);
local_properties = [cc,...                  1
                    gen',...                2
                    vul,...                 3
                    mean_vul_prey,...       4
                    mean_gen_pred,...       5
                    frac_basal_con,...      6
                    sTL,...                 7
                    paTL,...                8
                    betweenness,...         9
                    inLoop,...              10
                    pageRank...             11
                    ...
                    ];              

end