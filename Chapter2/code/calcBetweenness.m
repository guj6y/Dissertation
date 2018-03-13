function [cEB] = calcBetweenness(res,con)
%This function calculates the ecological betweenness for all consumers
%using Brande's algorithm.

n = max([res;con]);
sList = 1:n;
cEB = zeros(n,1);
canns = res==con;
res(canns) = [];
con(canns) = [];

for s = sList
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
        %We are at a new node if dist==inf
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



