function [avg] = wmean(x,w)

%wmean computes the average of x using weights w; by default this will omit
%nans in the thing being averaged.
x(w==0) = 0;
w(isnan(x)) = 0;
avg = sum(x.*w,'omitnan')./sum(w);

end