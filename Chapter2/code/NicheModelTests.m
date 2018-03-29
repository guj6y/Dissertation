%This constructs the six web models I will look at.

%first, need to get the empirical Data.
%Agglomeration level:
thisAggLevel = 41;
%Correspond to .41, the most agglomeration before we start to really lose
%accuracy in the classification tree.

load('selectors.mat');

%This is the average differences, std differences, etc. of all webs at each
%agglomeration level. Also includes the diffs for each web, just in case
%we want that, too.
T = load('groupedDistances.mat');
groupedDistanceMeans = T.groupedDistanceMeans;
groupedDistanceStd = T.groupedDistanceStd;
groupedDistanceN = T.groupedDistanceN;
groupedDistanceTs = T.groupedDistanceTs;
groupedDistancePs = T.groupedDistancePs;
allMeans= T.plotByMinDistance;

%We want the global properties only: need to get the Ss and Cs.
T = load('../../Chapter2/code/AgglomerationPropsMajority.mat','propsGlobal','globalCol');

propsGlobal = T.propsGlobal;
globalCol = T.globalCol;

clear T

