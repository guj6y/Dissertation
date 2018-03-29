function [connected] = walk(source,connected,res,con) 

%This function finds all nodes connected to the source using adjacency
%matrix 

neighbors = con(res==source);
connected = [connected;source];


if ~isempty(neighbors)
    neighbors(sum(neighbors'==connected)>0) = [];
end

if ~isempty(neighbors)
try
    for kk = neighbors'
        if sum(kk == connected)
            continue
        else
            connected = walk(kk,connected,res,con);
        end
    end
catch
    fpritnf('wtf')
end
end
   
end