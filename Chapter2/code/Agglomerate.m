function [propertiesMax,propertiesMean,propertiesMin] = Agglomerate(res,con,para0)



S = max([res;con]);

localPropertiesMax = nan(S,2,4);
localPropertiesMean = nan(S,2,4);
localPropertiesMin = nan(S,2,4);

globalPropertiesMax = nan(S,10);
globalPropertiesMean = nan(S,10);
globalPropertiesMin = nan(S,10);

globalPropertiesMax(:,1) = 1:S;
globalPropertiesMean(:,1) = 1:S;
globalPropertiesMin(:,1) = 1:S;

linkageMatrixCell = cell(S,1);
speciesPerNodeCell = cell(S,1);
parasitesPerNodeCell = cell(S,1);

maxAdjacencyMatrixCell = cell(S,1);
meanAdjacencyMatrixCell = cell(S,1);
minAdjacencyMatrixCell = cell(S,1);

S = max([res;con]);
res0 = res;
con0 = con;
linkageMatrix = full(sparse(res0,con0,1,S,S));


jaccardMatrices = cell(S-10,1);

Iall = zeros(S-1,1);
Jall = zeros(S-1,1);
%calculate the jaccard distance between all nodes:
jaccardMatrix = 1 - calculateSimilarity(res,con);

nSpeciesPerNode = ones(S,1);
nParasitesPerNode = para0+0;
ii = 0;
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
    
    %We preferentially merge smaller nodes first.
    [~, IDX] = sort(sourceSinkSize(:,1));
    sourceSinkSizeSorted = sourceSinkSize(IDX,:);
    IJSorted = IJSourceSinkSorted(IDX,:);
    
    %Now, preferentially merge pairs with same size source to the smaller
    %sink: tie breaker strives to equalize node size.
    for thisSourceSize = unique(sourceSinkSizeSorted(:,1))'
        %Find the links with the same source size:
        sortTheseLinks = sourceSinkSizeSorted(:,1) == thisSourceSize;
        
        %Pick out the subsets of size array and IJ array
        linksThisSourceSize = sourceSinkSizeSorted(sortTheseLinks,:);
        IJThisSourceSize = IJSorted(sortTheseLinks,:);
        
        %Figure out the right order for those subsets.
        [~, IDX] = sort(linksThisSourceSize(:,2));
        
        %Sort those subsets.
        linksThisSourceSizeSorted = linksThisSourceSize(IDX,:);
        IJThisSourceSizeSorted = IJThisSourceSize(IDX,:);
        
        %finally, need to break ties between sink sizes.
        for thisSinkSize = unique(linksThisSourceSizeSorted(:,2))'
            %Of the links with the given source size, pick those with sinks
            %of the current size.
            breakTheseSinkTies = linksThisSourceSizeSorted(:,2)==thisSinkSize;
            %How many links need to be re-arranged?
            numSinkTies = sum(breakTheseSinkTies);
            %Get a new order for those to be re-arranged.
            brokenTies = randsample(numSinkTies,numSinkTies);
            %rearrange those links only.
            IJThisSourceSizeSorted(breakTheseSinkTies,:) = IJThisSourceSizeSorted(brokenTies,:);
            linksThisSourceSizeSorted(breakTheseSinkTies,:) = linksThisSourceSizeSorted(brokenTies,:);
        end
        %Accept the new ordering of this portion of the IJsorted array.
        IJSorted(sortTheseLinks,:) = IJThisSourceSizeSorted;
        %likewise for thesourceSinkSizeSOrted array. (this might not be
        %used from her on out.
        sourceSinkSizeSorted(sortTheseLinks,:) = linksThisSourceSizeSorted;
        
    end
    
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
        speciesPerNodeCell{ii} = nSpeciesPerNode;
        parasitesPerNodeCell{ii} = nParasitesPerNode;
        linkageMatrixCell{ii} = linkageMatrix;
        jaccardMatrices{ii} = jaccardMatrix;
        
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

        nNodes = length(jaccardMatrix);
        if nNodes ==10
            break
        end
        
        nSpeciesPerNode(J(1)) = nSpeciesPerNode(J(1)) + nSpeciesPerNode(I(1));
        nParasitesPerNode(J(1)) = nParasitesPerNode(J(1)) + nParasitesPerNode(I(1));
        
        speciesPerNodeCell{ii} = nSpeciesPerNode;
        parasitesPerNodeCell{ii} = nParasitesPerNode;
        nSpeciesPerNode(I(1)) = [];
        nParasitesPerNode(I(1)) = [];
        
        para = nParasitesPerNode./nSpeciesPerNode;
        
        maxLinkageMatrix = nSpeciesPerNode*nSpeciesPerNode';
        maxA = linkageMatrix>0;
        meanA = linkageMatrix>=0.5*maxLinkageMatrix;
        minA = linkageMatrix==maxLinkageMatrix;
        
        
        
        maxAdjacencyMatrixCell{ii} = maxA;
        meanAdjacencyMatrixCell{ii} = meanA;
        minAdjacencyMatrixCell{ii} = minA;
        if max(size(maxA))<2
            fprintf('wtfmate?\n')
        end
        
        [localMaxProps, globalMaxProps] = calculateAgglomerateProperties(maxA,para);
        localPropertiesMax(ii,:,:) = localMaxProps;
        globalPropertiesMax(ii,2:end) = globalMaxProps;
        
        [localMeanProps, globalMeanProps] = calculateAgglomerateProperties(meanA,para);
        localPropertiesMean(ii,:,:) = localMeanProps;
        globalPropertiesMean(ii,2:end) = globalMeanProps;
        
        [localMinProps, globalMinProps] = calculateAgglomerateProperties(minA,para);
        localPropertiesMin(ii,:,:) = localMinProps;
        globalPropertiesMin(ii,2:end) = globalMinProps;


        I(I>I(1)) = I(I>I(1))-1;
        J(J>I(1)) = J(J>I(1))-1;
        I(1) = [];
        J(1) = [];

    end

end

propertiesMax.local = localPropertiesMax;
propertiesMax.global = globalPropertiesMax;

propertiesMean.local = localPropertiesMean;
propertiesMean.global = globalPropertiesMean;

propertiesMin.local = localPropertiesMin;
propertiesMin.global = globalPropertiesMin;

end
