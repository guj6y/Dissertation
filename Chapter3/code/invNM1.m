function [out] = invNM1(in,varargin)

    globalCol = varargin{1};
    S = in(globalCol('S'));
    Sf = in(globalCol('Sf'));
    Sp = S - Sf;
    C = in(globalCol('C'));
    [res,con,para] = inverseNicheModelC(S,C,Sp);
    
    out =  webPropsRaw(res,con,para);
end