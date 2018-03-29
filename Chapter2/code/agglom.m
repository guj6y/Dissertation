function [properties,fullStruct,reducedStruct] = agglom(res,con,para0,fLinkage)

S = numel(para0);

adjMatricesReduced = cell(S,1);
jacMatricesReduced = cell(S,1);

linkageMatrices = cell(S,1);
jacMatricesFull = cell(S,1);

paraFullCell = cell(S,1);
carnFullCell = cell(S,1);
freeFullCell = cell(S,1);

paraReducedCell = cell(S,1);
carnReducedCell = cell(S,1);
freeReducedCell = cell(S,1);

paraFullCell{1} = para0;
paraReducedCell{1} = para0;

nSpeciesPerNode = ones(S,1);
nParasitesPerNode = para0+0;
para0=para0+0;

A = sparse(res,con,1,S,S);
basal = full((sum(A)==0)');
carn0 = true(S,1);
carn0(con(basal(res))) = false;
carn0(basal) = false;
carn0(para0>0) = false;

linkageMatrices{1} = A;
carnFullCell{1} = carn0;
carnReducedCell{1} = carn0;

freeFullCell{1} = ~(basal|para0);
freeReducedCell{1} = ~(basal|para0);

carn0 = carn0+0;
nCarnivoresPerNode = carn0;
nBasalPerNode = basal;
[nodalProps,globalProps,localMeans] = agglomProps(A,para0,nSpeciesPerNode,0,1);
jaccardMatrix = 1 - calculateSimilarity(res,con);

[S,nLocalProps] = size(nodalProps);
nGlobalProps = length(globalProps);
nLocalMeans = length(localMeans);

bigLocalPropertiesMatrix = -9999*ones(S*(S+1)/2,nLocalProps);
bigGlobalPropertiesMatrix = -9999*ones(S,nGlobalProps);
bigLocalMeansMatrix = -9999*ones(S,nLocalMeans);


bigLocalPropertiesMatrix(1:S,1:nLocalProps) = nodalProps;
bigGlobalPropertiesMatrix(1,:) = globalProps;
bigLocalMeansMatrix(1,:) = localMeans;
adjMatricesReduced{1} = A;
jacMatricesFull{1} = jaccardMatrix;
jacMatricesReduced{1} = jaccardMatrix;
startingPoints = [0;cumsum((S:-1:2)')+1];
findStartingPoints = linspace(S,1,S)';

res0 = res;
con0 = con;
linkageMatrix = full(sparse(res0,con0,1,S,S));

Iall = zeros(S-1,1);
Jall = zeros(S-1,1);
%calculate the jaccard distance between all nodes:


ii = 1;
nNodes = length(jaccardMatrix);
while nNodes > 10
    minDistance = min(jaccardMatrix(:));
    [I, J] = find(triu(jaccardMatrix,1)==minDistance);
    
    if minDistance ==1
        break
    end
    IJ = [I J];
    [nRows,~] = size(IJ);
    %We think of merging the smaller node into the larger node; if they
    %are the same size we choose randomly. (smaller node = source,
    %larger node = sink).
    if nRows >1
        [sourceSinkSize,IDX] = sort(nSpeciesPerNode(IJ),2);
    else
        [sourceSinkSize,IDX] = sort(nSpeciesPerNode(IJ)',2);
    end
    ties = sourceSinkSize(:,1)==sourceSinkSize(:,2);
    nTies = sum(ties);
    tieBreakers = randi(2,nTies,1);
    IDX(ties,1) = tieBreakers;
    IDX(ties,2) = mod(tieBreakers,2)+1;
    
    IJSourceSinkSorted = IJ(nRows*(IDX-1)+(1:nRows)');
    sourceSinkSize = sourceSinkSize(nRows*(IDX-1)+(1:nRows)');
    tieBreakers = rand(size(I));
    %DID WAY TOO MUCH WORK HERE. Sortrows is way better!
    [~, Idx] = sortrows([sourceSinkSize tieBreakers], [2,1,3]);
    IJSorted = IJSourceSinkSorted(Idx,:);
    
    %NOW, re-define I and J:
    I = IJSorted(:,1);
    J = IJSorted(:,2);
   
    
    %Code that checks if we are okay here: Okay means only combine
    %"fully connected graphs"; strongly connected graphs that are not
    %fully connected are problematic, to say the least.
    %First question: are either of the first pair involved in any other
    %links?
    
    while numel(I)>0
        ii = ii+1;
        
        
        %The order of the links is now set. If the current source (I(1) =
        %Ihead) is supposed to merge with any other sinks, remove those
        %mergers.
        Ihead = I(1);
        Jhead = J(1);
        Itail = I(2:end);
        Jtail = J(2:end);
        
        Jtail(Itail==Ihead) = [];
        Itail(Itail==Ihead) = [];
        
        
        %IT is possible, In the case of source-sink size ties that a source
        %will later be a sink for a node of the same size. We want to
        %remove these ambiguities as well.
        Itail(Jtail==Ihead) = [];
        Jtail(Jtail==Ihead) = [];
        
        %Now, reform the I,J arrays:
        
        I = [Ihead;Itail];
        J = [Jhead;Jtail];
        
        %This is for looking at who is in each node; will be used in calculating
        %purity of nodes at each step; careful since this uses the updated
        %numbering at each step.
        Iall(ii) = I(1);
        Jall(ii) = J(1);

        %The adjacency matrix is shrinking; each entry i,j tells you how many links
        %exist inthe OG data from i to j. Note that there could be many links to a
        %single node in the data; I wanted to distinguish number of links from i
        %and the number of links to j
        likageMatrixNew = linkageMatrix;
        likageMatrixNew(:,J(1)) = linkageMatrix(:,J(1)) + linkageMatrix(:,I(1));
        likageMatrixNew(J(1),:) = linkageMatrix(J(1),:) + linkageMatrix(I(1),:);
        
        likageMatrixNew(I(1),:) = [];
        likageMatrixNew(:,I(1)) = [];

        linkageMatrix = likageMatrixNew;
        linkageMatrices{ii} = linkageMatrix;
        
        %Update the distance matrices
        jaccardMatrixNew = jaccardMatrix;

        newDistancesToJ1 = (jaccardMatrix(:,J(1))*nSpeciesPerNode(J(1)) + ...
            jaccardMatrix(:,I(1))*nSpeciesPerNode(I(1)))/(nSpeciesPerNode(J(1))+nSpeciesPerNode(I(1)));

        jaccardMatrixNew(:,J(1)) = newDistancesToJ1;
        jaccardMatrixNew(J(1),:) = newDistancesToJ1;
        jaccardMatrixNew(I(1),:) = [];
        jaccardMatrixNew(:,I(1)) = [];
        jaccardMatrixNew(eye(size(jaccardMatrixNew))>0)=1;
        jaccardMatrix = jaccardMatrixNew;
        
        jacMatricesFull{ii} = jaccardMatrix;
        
        
        nNodes = length(jaccardMatrix);
        
        nSpeciesPerNode(J(1)) = nSpeciesPerNode(J(1)) + nSpeciesPerNode(I(1));
        nParasitesPerNode(J(1)) = nParasitesPerNode(J(1)) + nParasitesPerNode(I(1));
        nCarnivoresPerNode(J(1)) = nCarnivoresPerNode(J(1)) + nCarnivoresPerNode(I(1));
        nBasalPerNode(J(1)) = nBasalPerNode(J(1)) + nBasalPerNode(I(1));
        
        nSpeciesPerNode(I(1)) = [];
        nParasitesPerNode(I(1)) = [];
        nCarnivoresPerNode(I(1)) = [];
        nBasalPerNode(I(1)) = [];
        
        %Fixing the counting after deleting nodes
        I(I>I(1)) = I(I>I(1))-1;
        J(J>I(1)) = J(J>I(1))-1;
        %Delete the nodes
        I(1) = [];
        J(1) = [];
        
        para = nParasitesPerNode./nSpeciesPerNode;
        carn = nCarnivoresPerNode./nSpeciesPerNode;
        basal = nBasalPerNode./nSpeciesPerNode;
        
        paraFullCell{ii} = para;
        carnFullCell{ii} = carn;
        freeFullCell{ii} = ~((basal>=0.5)|(para>=0.5));
        
        %Max Linkage Criterion
        %{
        A = linkageMatrix>0;
        adjacencyMatrixCell{ii} = A;
        %}
        maxLinkageMatrix = nSpeciesPerNode*nSpeciesPerNode';        
        %Mean linkage criterion
        %
        %
        A = linkageMatrix>=fLinkage*maxLinkageMatrix;
        
        %{
        %Min Linkage Criterion
        A = linkageMatrix==maxLinkageMatrix;
        %}
      
        
        [localProps, globalProps, localMeans, res, con] = agglomProps(A,para,nSpeciesPerNode,minDistance,ii);
        
        nNodesReduced = globalProps(1);
        startHere = startingPoints(nNodes==findStartingPoints);
        
        
        if res == 0 
            A = 0;
            jaccardReduced = 0;
            nNodesReduced = 1;
            endHere = startHere + nNodesReduced - 1;
        else
            A = sparse(res,con,1,nNodesReduced,nNodesReduced)>0;
            jaccardReduced = 1 - calculateSimilarity(res,con);
            endHere = startHere+nNodesReduced-1;
        end
        
        
        jacMatricesReduced{ii} = jaccardReduced;
        bigLocalPropertiesMatrix(startHere:endHere,:) = localProps;
        bigGlobalPropertiesMatrix(ii,:) = globalProps;
        bigLocalMeansMatrix(ii,:) = localMeans;
        adjMatricesReduced{ii} = A;
        
        paraReducedCell{ii} = localProps(:,10);
        carnReducedCell{ii} = localProps(:,11);
        freeReducedCell{ii} = ~((localProps(:,2)==0)|localProps(:,10));
        
        


    end

end


properties.local = bigLocalPropertiesMatrix(bigLocalPropertiesMatrix(:,1)~=-9999,:);
properties.global = bigGlobalPropertiesMatrix(bigGlobalPropertiesMatrix(:,1)~=-9999,:);
properties.localMeans = bigLocalMeansMatrix(bigLocalMeansMatrix(:,1)~=-9999,:);
jacMatricesFull(cellfun(@isempty,jacMatricesFull)) =[];
linkageMatrices(cellfun(@isempty,linkageMatrices)) =[];
paraFullCell(cellfun(@isempty,paraFullCell)) =[];
carnFullCell(cellfun(@isempty,carnFullCell)) =[];
freeFullCell(cellfun(@isempty,freeFullCell)) =[];
jacMatricesReduced(cellfun(@isempty,jacMatricesReduced)) =[];
adjMatricesReduced(cellfun(@isempty,adjMatricesReduced)) =[];
paraReducedCell(cellfun(@isempty,paraReducedCell)) =[];
carnReducedCell(cellfun(@isempty,carnReducedCell)) =[];
freeReducedCell(cellfun(@isempty,freeReducedCell)) =[];

fullStruct = struct('jac',{jacMatricesFull}...
                  ,'linkage',{linkageMatrices}...
                  ,'para',{paraFullCell}...
                  ,'carn',{carnFullCell}...
                  ,'free',{freeFullCell}...
                  );
reducedStruct = struct('jac',{jacMatricesReduced}...
                     ,'linkage',{adjMatricesReduced}... 
                     ,'para',{paraReducedCell}...
                     ,'carn',{carnReducedCell}...
                     ,'free',{freeReducedCell}...
                     );

%{
propertiesMean.local = localPropertiesMean;
propertiesMean.global = globalPropertiesMean;
propertiesMean.corr = corrMeanMatrices;

propertiesMin.local = localPropertiesMin;
propertiesMin.global = globalPropertiesMin;
propertiesMin.corr = corrMinMatrices;
%}


end
