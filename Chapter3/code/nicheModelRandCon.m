function [out] = nicheModelRandCon(in,varargin)
    globalCol = varargin{1};
    S = in(globalCol('S'));
    C = in(globalCol('C'));
    fPar = in(globalCol('fPar'));
    [res, con] = NicheModel_nk(S,C);
    
    basal = full(sum(sparse(res,con,1,S,S))==0)';
    SList = 1:S;
    paraIdx = randsample(SList(~basal),round(fPar*sum(~basal)));
    para = false(S,1);
    para(paraIdx) = true;
    
    out =  webPropsRaw(res,con,para);
end