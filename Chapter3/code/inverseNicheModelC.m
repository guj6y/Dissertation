function [res,con,para,n_new,c_new,r_new]= ...
     inverseNicheModelC(S, C, Sp)                
%------------------------------------------------------------------
Sf = S - Sp;
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
    tries=tries-1;
    web_connected = false;
    %assign niche values from a uniform distribution
    n = rand(S,1);
    para = false(S,1);
    %his introduces some minor bias to the parasites. should probably fix
    %this. (actually, easy to do! it is uniformly distributed between 0 and
    %the max (at least, approximately...))
    
    para(randsample(SList(n<max(n)),Sp)) = true;
    free = ~para;    
    Enp = max(n)/2;
    %designate range for each species

    %parameters for beta distribution:
    alpha = 1;
    
    %Parasites will always have a host, so we have to correct their C
    %accordingly.
    Cpara = (C*S^2 - C*Sf*S - Sp)/((S-1)*Sp);
    Cfree = (C*S^2 - C*Sp*S)/((Sf)*S);
    betaPara = (1-Enp)/Cpara - 1;
    betaFree = (1-2*Cfree)/(2*Cfree); 
    
    r = zeros(S,1);
    r(free) = betarnd(alpha,betaFree,Sf,1); 
    r(para) = betarnd(alpha,betaPara,Sp,1); 
    
    %vector of ranges: 
    r(free) = r(free).*n(free);
    r(para) = r(para).*(1-n(para));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set center of range, uniformly distributed in [r_i/2,n_i];
    c = zeros(S,1);
    c(free) = min(1-r(free)./2,rand(Sf,1).*(n(free)-r(free)./2)+r(free)./2);
    c(para) = cellfun(@(x) randsample(n(n>x),1),num2cell(n(para)));
    c((c + r/2)>1) = 1-r(c+r/2>1)/2;
    c((c - r/2)<0) = r(c-r/2<0)/2;

    r(n==min(n(free))) = 0;
    %change the r of highest index species to 0
    %so we have a basal species in every web
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    c_new = c;
    r_new = r;
    n_new = n;
    
    
    %lower border of niche range for every prey:
    preymins = c_new - r_new/2;
    
    %upper border of niche range for every predator:
    preymaxs = c_new + r_new/2;

    %fills the empty matrix with niche ranges:
    n_mx = n_new*ones(1,S); 
    %matrix with the lowest points of ranges in every column:
    preymins_mx = ones(S,1)*preymins'; 
    %same, with highest:
    preymaxs_mx = ones(S,1)*preymaxs';

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
    
    %'Walk' is where most of the time is spent.  At this point, however, it
    %takes about .04s to run on a 100 species, .12 connectance web.  IF you
    %need to, you could use the shortest paths matrix from floyd warshall
    %to figure out if everything goes to a basal.  Could maybe also figure
    %out if it is connected using the sp stuff - one function call for both
    %things might save some time.
    
        % indices of links (which element in web_mx?)
    links = numel(res);  
    % Actual connectance
    C_web = links/(S^2);  
    if (abs(C_web-C)*1.0/C) > error_n
        ok_n=0;
    else
        ok_n=1;
    end

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

    
end