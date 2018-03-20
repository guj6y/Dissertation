function [pr] = pageRank2(res,con)

%This function calculates an ecological pagerank centrality meaure
%(Allesina2009). This is called an eigenvector based measure because the
%values for each node are the entries of the dominant eigenvector for the
%new matrix. Since we don't want to return the value of the root node, we
%remove it from the entries and then normalize so that the sum is equal to
%one.


canns = res==con;
res(canns) = [];
con(canns) = [];
n = max([res;con]);
sList =1:n;

A = sparse(res,con,1,n,n);
basal = (sum(A)==0)';
basalIds = sList(basal);

resNew = [(1:n)';repmat(n+1,sum(basal),1)];
conNew = [repmat(n+1,n,1);basalIds'];

res_ = [res;resNew];
con_ = [con;conNew];

A = sparse(res_,con_,1,n+1,n+1);
S = A./sum(A,2);
try
[pr,~] = eigs(S',1);
catch
    pr = nan(n,1);
end
pr = abs(pr(1:n))/norm(pr(1:n),1);



end