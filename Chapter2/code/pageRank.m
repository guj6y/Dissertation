function [pr] = pageRank(res,con)

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

A = full(sparse(res_,con_,1,n+1,n+1));
%OG Allesina Paper:
S = A./sum(A);

%Naive Markov Chain idea (nutrients from resources evenly split FROM the 
%resource node).
%S2 = (A./sum(A,2))';

%Better Markov Chain Idea (nutrients from resources split according to
%generality of the consumers. These are transposed so that the proper
%directionality is espected (right vs. left eigenvectors).
%S = A./sum(A);
%S = (S./sum(S,2))'; 

%Use perron-frobenius theorem to find the eigenvector to dominant
%eigenvalue of one; what makes this nice.
%try
    [pr,~] = eigs(S,1);
    pr = abs(pr(1:end-1))/norm(pr(1:end-1),1);
%catch
%    pr = nan(n,1);
%end
    




end