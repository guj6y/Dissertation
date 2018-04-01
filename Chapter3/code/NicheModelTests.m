function NicheModelTests(distIndex,linkageType)
%This constructs the six web models I will look at.
%Correspond to .41, the most agglomeration before we start to really lose
%accuracy in the classification tree.

%distancesPlotted = [0 0.05 0.1 0.2 0.4];
addpath(genpath('~/matlab_bgl'));

minDistances = [0,0.05,0.1,0.2,0.4];
minDistance = minDistances(distIndex);
nWebsPerWeb = 1000;

webFunctions = {@nicheModelRandCon;
        ...@nicheModelRandCarn;
        ...@nicheModelRandTree;
        @invNM1;
        @invNM2;
        @invNM3};


nModels = numel(webFunctions);
nWebs = 6;
nLocalProps = 12;
nGlobalProps = 18;

%for linkages = {'Min','Max'}

%linkageType = linkages{:};

load(sprintf('../../Chapter2/code/selectors%sLinkage.mat',linkageType));


T = load(sprintf('../../Chapter2/code/AgglomerationProps%sLinkage.mat',linkageType));

propsGlobal = T.propsGlobal;
globalCol = T.globalCol;
reducedProps = T.reducedProps;


clear T

adjMatrices = {reducedProps{1}.linkage{:}...
          ,reducedProps{2}.linkage{:}...
          ,reducedProps{3}.linkage{:}...
          ,reducedProps{4}.linkage{:}...
          ,reducedProps{5}.linkage{:}...
          ,reducedProps{6}.linkage{:}...
         }';



localVariableNames = {'$v$'...          1
                     ,'$g$'...          2
                     ,'$T$'...         3    
                     ,'$Fr$'...           4
                     ,'$C_B$'...         5
                     ,'$C_{EB}$'...     6
                     ,'$v_r$'...  7
                     ,'$g_c$'...  8
                      ,'$\gamma^{cyc}$'...    14
                      ,'$\gamma^{mid}$'...    15
                      ,'$\gamma^{in}$'...     16
                      ,'$\gamma^{out}$'...    17
                     };


count = 0;
%for minDistance = distancesPlotted
count = count+1;
%first, need to get the empirical Data.
%Agglomeration level:

getEm = find(diff(propsGlobal(:,globalCol('d_J'))<=minDistance)==-1);
empPara = cellfun(@(x) x.para(:)',reducedProps,'UniformOutput',false)';
empPara = [empPara{:}]';
empPara = empPara(getEm);
empCarn = cellfun(@(x) x.carn(:)',reducedProps,'UniformOutput',false)';
empCarn = [empCarn{:}]';
empCarn = empCarn(getEm);

empParaCell = empPara;
empPara = cell2mat(empParaCell);

empWebs = adjMatrices(getEm);

resCell = cell(6,1);
conCell = cell(6,1);
for ii = 1:6
    [resCell{ii},conCell{ii}] = ind2sub(size(empWebs{ii}),find(empWebs{ii}));
end

propertiesWebsToMatch = propsGlobal(getEm,:);


%The plan is to calculate 100 webs with each proposed model for each web.
%We then compare the average performance across all 6 webs.

globalProps = cell(nModels,nWebs);
freeProps = cell(nModels,nWebs);
paraProps = cell(nModels,nWebs);

[globalProps{:}] = deal(zeros(nWebsPerWeb,nGlobalProps));
[freeProps{:}] = deal(zeros(nWebsPerWeb,nLocalProps));
[paraProps{:}] = deal(zeros(nWebsPerWeb,nLocalProps));

%Get the empirical properties; after the next two lines,
%empiricalproperties is 2x6; first row is global properties, second row is
%local.
empiricalProperties = cellfun(@(x,y,z) webPropsRaw(x,y,z>0),resCell,conCell,empParaCell,'UniformOutput',false);
empiricalProperties = reshape([empiricalProperties{:}],3,6);

nCores = feature('numcores');
parpool('local',nCores/2);

for ii = 1:nModels

   webFunction = webFunctions{ii};
   fprintf('nicheModel no %u\n',ii);
   
   for jj = 1:nWebs
      
       %Just takes it all in; wanted this to be as standardized as
       %possible :):
       fprintf('web no %u\n',jj);
       webFunctionIn = propertiesWebsToMatch(ii,:);
       tempGlobal = zeros(nWebsPerWeb,nGlobalProps);
       tempFree = zeros(nWebsPerWeb,nLocalProps);
       tempPara = zeros(nWebsPerWeb,nLocalProps);
      
       parfor kk = 1:nWebsPerWeb
            
           %WebFunctionOut is a cell array with the appropriate carn-para
           %differences as well as the global properites.
           webFunctionOut = webFunction(webFunctionIn,globalCol);

           tempGlobal(kk,:) = webFunctionOut{1};
           tempFree(kk,:) = webFunctionOut{2};
           tempPara(kk,:) = webFunctionOut{3};
       
       end
       globalProps{ii,jj} = tempGlobal;
       freeProps{ii,jj} = tempGlobal;
       paraProps{ii,jj} = tempLocal;
       warning('off','all');
   
   end

end

fprintf('saving...\n');
save(sprintf('NicheTestResults%sLinkage-Distance%u',linkageType,distIndex)...
    ,'globalProps','freeProps','paraProps');
fprintf('done.\n');
clear carnParaDiffs globalProps treeProps
save(sprintf('nicheModelTestsVariables%Linkage-Distance%u.mat',linkageType,distIndex));

end
