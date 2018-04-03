function [cEB,spsFromBasal,totalDistanceFromBasal,minSpToBasal,numBasalCon] = ecoBtwn(res,con)
%This function calculates the ecological betweenness for all consumers
%using Brande's algorithm.

n = max([res;con]);
sList = 1:n;
cEB = zeros(n,1);
minSpToBasal = inf(n,1);
numBasalCon = zeros(n,1);
spsFromBasal = zeros(n,1);
totalDistanceFromBasal = zeros(n,1);
basalBin = false(n,1);

for ii = 1:n
    if sum(ii == con)>0
        continue
    else
        basalBin(ii) = true;
    end
    
end
basalID = find(basalBin)';
nBasal = numel(basalID);
canns = res==con;
res(canns) = [];
con(canns) = [];
coder.varsize('S','Q','Qnew','Snew',10000);

for ii = 1:nBasal
    s = basalID(ii);
    S = zeros(0,1);
    pred = zeros(n,n);
    sig = zeros(n,1);
    sig(s) = 1;
    dist = inf(n,1);
    dist(s) = 0;
    Q = s;
    
    while numel(Q)>0
        v = Q(1);
        a = numel(Q);
        Qnew = zeros(a-1,1);
        Qnew(:) = Q(2:end);
        S = [v; S];
        w = con(res==v);
        %We are at anew node if dist==inf
        newNodes = isinf(dist(w));
        dist(w(newNodes)) = dist(v)+1;
        nQAdd = w(newNodes);
        Q = [Qnew;w(newNodes)];


        %We found a new sp to w if dist(w) == dist(v) +1.
        vPredOf = w(dist(w)==(dist(v)+1));
        numPreds = numel(vPredOf);
        
        
        if numPreds>0
            changeMe = pred(vPredOf,:);
            changeMe((sum(changeMe>0,2))*numPreds + (1:numPreds)') = v;
            pred(vPredOf,:) = changeMe;
            sig(vPredOf) = sig(vPredOf) + sig(v);
        end
    end
    spsFromBasal = spsFromBasal + sig;
    pathsExist = isfinite(dist);
    totalDistanceFromBasal(pathsExist) = totalDistanceFromBasal(pathsExist) + dist(pathsExist).*sig(pathsExist);
    minSpToBasal = min(minSpToBasal,dist);
    numBasalCon = numBasalCon + 1*(isfinite(dist)&(dist>0));
    
    delta = zeros(n,1);
    
    while numel(S) > 0
        w = S(1);
        a = numel(S);
        Snew = zeros(a-1,1);
        Snew(:) = S(2:end);
        S = Snew;
        pred_w = pred(w,pred(w,:)>0);
        nPred_w = numel(pred_w);
        for jj = 1:nPred_w
            v = pred_w(jj);
            delta(v) = delta(v) + (sig(v)/sig(w))*(1+delta(w));
        end
        if w~=s
            cEB(w) = cEB(w) + delta(w);
        end
    end
        
    
end

end


