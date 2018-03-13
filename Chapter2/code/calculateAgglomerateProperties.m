function [localProperties,globalProperties] = calculateAgglomerateProperties(A,para)

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
    if sizeii>maxCompSize
        maxCompSize=sizeii;
        largeComp = ii;
    end
end
numLargest = sum(compSizes==maxCompSize);

if numLargest>1
    if numLargest==S0
        fprintf('Completely disconnected. Returning nans.\n')
        localProperties = nan(2,4);
        globalProperties = nan(1,9);
        return
    end
    comps = 1:n;
    fprintf('Tie in largest component... averaging them.\n')
    localProperties = zeros(2,4);
    globalProperties = zeros(1,9);

    for ii = comps(compSizes==maxCompSize)
        B = A(grps==ii,grps==ii);
        para_ = para(grps==ii);
        [loc,glob] = calculateAgglomerateProperties(B,para_);
        localProperties = localProperties+loc/numLargest;
        globalProperties = globalProperties+glob/numLargest;
    end
    globalProperties(end) = n;
    return
else
    %I like winners. Winners get to stay.
    winners = grps==largeComp;
    A = A(winners,winners);
    para = para(winners);
    
    %Now, we calculate properties using the largest connected component.
    S = length(A);
    L = sum(sum(A));

    %In degree (generality)
    gen = full(sum(A))/(L/S);
    gen = gen';

    %This defines parasite cutoffs; without this line we calculate weighted
    %averages:
    %para = para>0.5;

    basal = gen==0;
    free = 1-para;
    free(basal) = 0;
    con = (~basal);

    genCon = sum(gen.*con)/sum(con);

    genPara = sum(gen.*para/sum(para))./genCon;
    genFree = sum(gen.*free/sum(free))./genCon;

    %out degree (vulnerability)
    vul = full(sum(A,2))/(L/S);
    vulCon = sum(vul.*con)./sum(con);

    vulPara = sum(vul.*para/sum(para))./vulCon;
    vulFree = sum(vul.*free/sum(free))./vulCon;

    %Formula for the prey-averaged trophic level.
    A(((1:S)-1)*S+(1:S)) = 0; %no cannibalism in this matrix; think it's better this way.
    patl_mx = sparse(A)*(diag(1./sum(sparse(A))));

    patl = (speye(S)-patl_mx')\ones(S,1);
    patlCon = sum(patl.*con)/sum(con);

    patlPara = sum(patl.*para/sum(para))./patlCon;
    patlFree = sum(patl.*free/sum(free))./patlCon;

    %Clustering Coefficient(s)
    try
    cc = clustering_coefficients(sparse(A));
    ccCon = sum(cc.*con)./mean(cc(con));
    catch
        fprintf('wtf, mate?\n')
    end
    ccPara = sum(cc.*para/sum(para))./ccCon;
    ccFree = sum(cc.*free/sum(free))./ccCon;

    localProperties = [vulFree, genFree, patlFree, ccFree;
                       vulPara, genPara, patlPara, ccPara];

    Lcon = sum(sum(A(con,con)));
    Scon = sum(con);
    Ccon = Lcon/Scon^2;

    para = para>=0.5;
    free = free>0.5;

    Lff = sum(sum(A(free,free)));
    Sf = sum(free);
    Cff = Lff/Sf^2;

    Lpp = sum(sum(A(para,para)));
    Sp = sum(para);
    Cpp = Lpp/Sp^2;

    Lfp = sum(sum(A(free,para)));
    Cfp = Lfp/(Sp*Sf);

    Lpf = sum(sum(A(para,free)));
    Cpf = Lpf/(Sp*Sf);

    fPar = Sp/(Sp+Sf);

    globalProperties = [S,Ccon,Cff,Cpp,Cfp,Cpf,fPar,numLargest,n];

end

end