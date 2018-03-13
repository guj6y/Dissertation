function [avg] = wmean(x,w)

%wmean computes the average of x using weights w.
x(w==0) = 0;
avg = sum(x.*w)./sum(w);

end