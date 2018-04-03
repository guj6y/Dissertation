function [out] = nicheModelRandTree(in,varargin)
    out = cell(2,1);
    globalCol = varargin{1};
    tree = varargin{2};
    
    S = in(globalCol('S'));
    C = in(globalCol('C'));
    fPar = in(globalCol('fPar'));
    [res, con] = NicheModel_nk(S,C);
    
    basal = full(sum(sparse(res,con,1,S,S))==0)';
    carn = true(S,1);
    carn(con(basal(res))) = false;
    carn(basal) = false;
    
    nPara = round(fPar*sum(~basal));
    X = propsForClassTree(res,con);
   
    carn = ~basal;
    nCarn = sum(carn);
    carnIdx = find(carn);
    X = X(carn,:);

    nParaSoFar = 0;
    carnAsPara = false(nCarn,1);
    [~,probPara] = tree.predict(X);
    probPara = probPara(:,2);
   while nParaSoFar<nPara
        draw = rand(nCarn,1);
        IDdPara = draw <= probPara;
        carnAsPara(IDdPara) = true;
        nParaSoFar = sum(carnAsPara);

        %Need to 'take care' of extra parasites if we have them.
        if nParaSoFar > nPara
            nToKill = nParaSoFar-nPara;
            [~,Idx] = sort(probPara);
            notPara = find(IDdPara(Idx),nToKill);
            carnAsPara(Idx(notPara))=false;
        end
        
        probPara(IDdPara) = 1;
   end

    para = false(S,1);
    para(carnIdx(carnAsPara)) = true;
    
    out =  webPropsRaw(res,con,para);
end