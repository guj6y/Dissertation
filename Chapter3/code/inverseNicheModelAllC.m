function [res,con,para,n_new,c_new,rf,rp]= ...
     inverseNicheModelAllC(S, Sp, C, Cff, Cpp, Cpf, Cfp)                
%------------------------------------------------------------------
Sf = S - Sp;
Cs = [];
Lff = Cff*Sf^2;
Lpp = Cpp*Sp^2;
Lpf = Cpf*Sp*Sf;
Lfp = Cfp*Sp*Sf;
L = Lff + Lpp + Lpf + Lfp;

%globalStream = RandStream.getGlobalStream;
%reset(globalStream);
SList = 1:S;
tries=10000;                                                               
%error on the connectances
error_n=0.03;

%validation of the plain niche web:
ok_n=0; 
web_connected = false;
while ((tries>0) && ~(ok_n&&web_connected))
    web_connected = false;
    ok_n = 0;
    tries=tries-1;
    %assign niche values from a uniform distribution
    n = rand(S,1);
    para = false(S,1);
    
    %his introduces some minor bias to the parasites. should probably fix
    %this. (actually, easy to do! it is uniformly distributed between 0 and
    %the max (at least, approximately...))
    %Does it affect free-livers.. it's randomly sampling uniform random
    %samples... 
    para(randsample(SList(n<max(n)),Sp)) = true;
    free = ~para;    
    Enp = max(n)/2;
    %designate range for each species

    %parameters for beta distribution:
    alpha = 1;
    %Parasites will always have a host, so we have to correct their C
    %accordingly. (Number of links is due to being parasite (one per
    %parasite) and stochastic (the rest).
    Cfp = (L - Lff - Lpf - Lpp - Sp)/((Sf-1)*Sp);
    
    betaFF = (1-2*Cff)/(2*Cff); 
    betaPP = (1-Enp)/Cpp - 1; 
    betaPF = (1-2*Cpf)/(2*Cpf);
    betaFP = (1-Enp)/Cfp - 1;
    
    %rp is the ranges for preying on parasites; rf is for free-livers. 
    %Could also hvae different centers, in principle, but I'm not going 
    %there.
    
    rp = zeros(S,1);
    rf = zeros(S,1);
    
    rp(free) = betarnd(alpha,betaPF,Sf,1); 
    rp(para) = betarnd(alpha,betaPP,Sp,1); 
    
    rf(free) = betarnd(alpha,betaFF,Sf,1); 
    rf(para) = betarnd(alpha,betaFP,Sp,1); 
    
    %vector of ranges: 
    rf(free) = rf(free).*n(free);
    rf(para) = rf(para).*(1-n(para));
    
    rp(free) = rp(free).*n(free);
    rp(para) = rp(para).*(1-n(para));
    
    rMax = max([rf,rp],[],2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set center of range, uniformly distributed in [r_i/2,n_i];
    c = zeros(S,1);
    c(free) = min(1-rMax(free)./2,rand(Sf,1).*(n(free)-rMax(free)./2)+rMax(free)./2);
    c(para) = cellfun(@(x) randsample(n(n>x),1),num2cell(n(para)));
    c((c + rMax/2)>1) = 1-rMax(c+rMax/2>1)/2;
    c((c - rMax/2)<0) = rMax(c-rMax/2<0)/2;

    %correct for this bias: you lose out on one free liver by default.
    rf(n==min(n(free))) = 0;
    rp(n==min(n(free))) = 0;
    %change the r of highest index species to 0
    %so we have a basal species in every web
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    c_new = c;
    n_new = n;
    
    %lower border of niche range for every prey:
    preyMinsF = c_new - rf/2;
    preyMinsP = c_new - rp/2;
    
    %upper border of niche range for every predator:
    preyMaxsF = c_new + rf/2;
    preyMaxsP = c_new + rp/2;

    %fills the empty matrix with niche ranges:
    n_mx = n_new*ones(1,S); 
    %matrix with the lowest points of ranges in every column:
    preymins_mx = free*preyMinsF' + para*preyMinsP'; 
    %same, with highest:
    preymaxs_mx = free*preyMaxsF' + para*preyMaxsP';

    %Construct the web adjacency matrix;
    %if species in the row is in the diet range of the species in the 
    %column, it gets eaten (if species i is eaten by species j, set (i,j) =
    %1).
    
    web_mx=(n_mx>=preymins_mx)&(n_mx<=preymaxs_mx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  if there is an isolated species                                    %
    %  or something not connected to a basal                               %
    %  or a disconnected group                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    

    
    [res,con] = find(web_mx);
    %make sure that all species have a directed connection from a basal
    %species; only if connected.
    
    % indices of links (which element in web_mx?)
    links = numel(res);  
    % Actual connectance
    C_web = links/(S^2);  
    C_fp = sum(free(res)&para(con))/(Sf*Sp);
    C_ff = sum(free(res)&free(con))/Sf^2;
    C_pf = sum(para(res)&free(con))/(Sf*Sp);
    C_pp = sum(para(res)&para(con))/Sp^2;
    if (abs(C_web-C)*1.0/C) > error_n
        ok_n=0;
        Cs = [Cs;C_fp];
    else
        ok_n=1;
    end


    
    %'Walk' is where most of the time is spent.  At this point, however, it
    %takes about .04s to run on a 100 species, .12 connectance web.  IF you
    %need to, you could use the shortest paths matrix from floyd warshall
    %to figure out if everything goes to a basal.  Could maybe also figure
    %out if it is connected using the sp stuff - one function call for both
    %things might save some time.
    
    if ok_n
        
        %Make sure the graph is weakly connected (no isolated species)
        weak_comp = graphconncomp(sparse(web_mx),'Directed',1,'Weak',1);
        web_connected = (weak_comp == 1);
        
        if ~web_connected
            continue
        end
        
        unconnected_species = (1:S)';
        basal_species = unconnected_species(sum(web_mx)==0)';
        connected_species = [];

        for kk = basal_species
            connected_species = walk(kk,connected_species,res,con);
        end

        unconnected_species(connected_species) = [];
        
        if isempty(unconnected_species)
            web_connected = 1;
        else
            web_connected = 0;
        end
    end
    
    
    
    

end

    [res,con] = find(web_mx);

    
end