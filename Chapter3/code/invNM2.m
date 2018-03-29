function [out] = invNM2(in,varargin)

    globalCol = varargin{1};
    S = in(globalCol('S'));
    Sf = in(globalCol('Sf'));
    Sp = S - Sf;
    C = in(globalCol('C'));
    Cf = in(globalCol('Cf'));
    [res,con,para] = inverseNicheModelCfCp(S,Sp,C,Cf);
    
    out =  webPropsRaw(res,con,para);
end