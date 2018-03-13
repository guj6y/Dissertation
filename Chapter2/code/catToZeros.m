function [x] = catToZeros(x,y)
%x_ = x;
    [height,~] = size(x);
    %x([zeros(height,1) diff((x==0)')']>0) = y;
    %x(diff(([ones(height,1) x]==0)')'>0) = y;
    x((sum(x>0,2))*height + (1:height)) = y;
end