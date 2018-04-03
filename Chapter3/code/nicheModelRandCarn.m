function [out] = nicheModelRandCarn(in,varargin)
    globalCol = varargin{1};
    S = in(globalCol('S'));
    C = in(globalCol('C'));
    fPar = in(globalCol('fPar'));
    [res, con] = NicheModel_nk(S,C);
    
    basal = full(sum(sparse(res,con,1,S,S))==0)';
    carn = true(S,1);
    carn(con(basal(res))) = false;
    carn(basal) = false;
    nCarn = sum(carn);
        
    SList = 1:S;
    nPara = round(fPar*sum(~basal));
    
    if nCarn >= nPara
        paraIdx = randsample(SList(carn),nPara);
    else
        paraIdx = randsample(SList(~basal),nPara);
        fprintf('Not enough carnivores\n');
    end
    
    para = false(S,1);
    para(paraIdx) = true;
    
    out =  webPropsRaw(res,con,para);
    
end