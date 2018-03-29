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

basalID = sList(sum(sparse(res,con,1,n,n))==0)';

canns = res==con;
res(canns) = [];
con(canns) = [];

for s = basalID'
    S = [];
    pred = zeros(n,n);
    sig = zeros(n,1);
    sig(s) = 1;
    dist = inf(n,1);
    dist(s) = 0;
    Q = s;
    
    while numel(Q)>0
        v = Q(1);
        Q(1) = [];
        S = [v; S];
        w = con(res==v);
        %We are at anew node if dist==inf
        newNodes = isinf(dist(w));
        dist(w(newNodes)) = dist(v)+1;
        Q = [Q; w(newNodes)];
        
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
        S(1) = [];
        pred_w = pred(w,pred(w,:)>0);
        for v = pred_w
            delta(v) = delta(v) + (sig(v)/sig(w))*(1+delta(w));
        end
        if w~=s
            cEB(w) = cEB(w) + delta(w);
        end
    end
        
    
end

end

