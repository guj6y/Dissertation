function [out] = invNM3(in,varargin)

    globalCol = varargin{1};
    S = in(globalCol('S'));
    Sf = in(globalCol('Sf'));
    Sp = S - Sf;
    C = in(globalCol('C'));
    Cff = in(globalCol('Cff'));
    Cpp = in(globalCol('Cpp'));
    Cpf = in(globalCol('Cpf'));
    Cfp = in(globalCol('Cfp'));
    
    [res,con,para] = inverseNicheModelAllC(S,Sp,C,Cff,Cpp,Cpf,Cfp);
    
    out =  webPropsRaw(res,con,para);
end