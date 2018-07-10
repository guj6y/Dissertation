%This code sets up parameters structures for a parallel implementation of
%the parasites experiment, then calls the actual integrating function and
% also parses the data.

%Setting up the arrays and structures.
S=40;

try 
    load '../raw/metaSimData.mat'
    load '../raw/simParams.mat'
    load '../raw/webData.mat'
 
    appending = true;
    fprintf('We are appending data to a simulation that has already been run.\n');

catch     
    appending = false;
    fprintf('We are doing a brand-new simulation!\n')
end 

%Code can handle appending new values of factors or starting fresh.
%This sets up the structure of the experiments themselves. 
if appending 
    kFreesDone = kFrees;
    nKFreesDone = 1:numel(kFrees);
    
    kParasDone = kParas;
    nKParasDone = 1:numel(kParas);
    
    nWebsDone = nWeb;
    
    fParDone = fParAll0;
    nFParDone = numel(fParAll0);
    
    %Allow for testing new free living and parasitic body size ratios. This is complicated enough;
    %I could imagine ways to allow for new dynamical models, new factors, new fractions of paras-
    %ites, more webs... but I don't think it will come to that, and Don't want to take the time to
    %figure that out.
    %All this could mix up the out.mat arrays, BUT they get processed in that linear order anyway.
    
    newKFrees = 3;
    newKParas = [-5];
    fParNew = [];
    nFParNew = numel(fParNew);
    nWebNew = 0;

    kFrees = [kFreesDone,newKFrees];
    kParas = [kParasDone,newKParas];
    
    nFact1 = numel(kFrees);
    nFact2 = numel(kParas);
    nFact3 = 2;
    nFact4 = 2;
    nFacts = [nFact1 nFact2 nFact3 nFact4];
    
    allModels = fullfact(nFacts);
    nAllModels = length(allModels);
    
    doneModels = sum(allModels(:,1)==nKFreesDone,2)&sum(allModels(:,2)==nKParasDone,2);
    modelsToRun = allModels(~doneModels,:);
    nNewModels = length(modelsToRun); 
    fParAll0 = [fParDone, fParNew];
    nFParAll = numel(fParAll0);
    nWeb = nWebsDone + nWebNew;

    nSims = nWeb*(nFact1-numel(nKFreesDone))... Run the new kFrees on all websat 0 percent parasites.
        + nWeb*nAllModels*nFParNew ... All models and all webs need to be run on the new parasite fractions
        + nWebNew*nAllModels*(nFParDone-1) ... All models and old parasite fractions run on the new webs.
        + nWebsDone*nNewModels*(nFParDone-1) ... New Models on Old Webs with Old Fractions.
        + nWebNew*(numel(nKFreesDone)); ... Old kFrees on new webs.

        newWebData = cell(nWebNew,1);
[newWebData{:}] = deal(struct('web',0 ...
                    ,'LL',zeros(0,2) ...
                    ,'B0',zeros(S,1) ...
                    ,'gr',zeros(S,1) ...
                    ,'patl',zeros(S,1) ...
                    ,'parOrder',zeros(S,1) ...
                    ));
webData = [webData; newWebData];
else %If we're not appending
    nWeb = 100;
    
    nBasal = 6;
    nFree = S-nBasal;
    
    kFrees = [1 2 3];    %Models(1): BSR exponents for free livers
    kParas = [-3 -4 -5];  %Models(2): BSR exponents for free livers
    fracFrees = [false true]; %Models(3): including fraction of free living time (binary)
    fracParas = [false true]; %Models(4): including concomittant links (binary)

    nFact1 = numel(kFrees);
    nFact2 = numel(kParas);
    nFact3 = numel(fracFrees);
    nFact4 = numel(fracParas);

    nFacts = [nFact1,nFact2,nFact3,nFact4];

    allModels = fullfact(nFacts);
    nModels = length(allModels);

    fParAll0 = [0 linspace(1,17,9)/nFree];
    nfPar = numel(fParAll0);

    nSims = nWeb*(nFact1) ... Run 1 simulation for each kFree @ 0% parasites.
          + nWeb*nModels*(nfPar-1); ... run each web for each model for all but one fraction of parasites.
    
    nWebNew = 0;
    nWebsDone = 0;
    nKFreesDone = [];
    kFreesDone = [];
    nFParDone = 0;
    nKParasDone = [];
    kParasDone = [];
end 

save('../raw/metaSimData.mat','S','nWeb','kFrees','kParas','fracFrees','fracParas','fParAll0','nFacts','nBasal')

simParams = cell(nSims,1);
[simParams{:}] = deal(struct('web',0 ...
    ,'fPar',0 ...
    ,'kFree',0 ...
    ,'kPara',0 ...
    ,'fracFree',0 ...
    ,'fracPara',0 ...
    ,'modelCode',[0 0]...
    ,'para',zeros(40,0)...
    ,'LL',[]...
    ,'B0',zeros(40,1)...
    ,'gr',zeros(40,1)...
    ,'parOrder',zeros(40,1)...
    ,'patl',zeros(40,1)...
    ));

%%% Parameters of the Niche Webs
S = 40;
C = .15;

simNo = 0;

%This populates the necessary web data for each simulation. The parameters that vary according to the
%needs of the factors.

for ii = 1:nWeb 
    %Generate nichewebs or load old one; save the web data separately for ease of access. Will also
    %get saved in the structure used by theparfor loop. Slight redundancy, but the memory usage pales
    %comparison to the time series data, so I think it is okay. It makes life much easier later on. 
    if ii <= nWebsDone
        B0 = webData{ii}.B0;
        res = webData{ii}.LL(:,1);
        con = webData{ii}.LL(:,2);
        gr = webData{ii}.gr;
        patl = webData{ii}.patl;
        idxPar = webData{ii}.parOrder;
        basal = gr>0;
        nFree = S-sum(basal);
    else
        webBad = true;
 
        while webBad
            [res, con,~,~,~] = NicheModel_nk(S,C);
            simMx = calculateSimilarity(res,con);
            mx = sparse(res,con,1,S,S);
            
            %Need to decide if each species is a basal.
            basal = false(S,1);
            
            for kk = 1:S
                if sum(con==kk)==0
                    basal(kk) = true;
                end
            end
            
            gr = basal;
            webBad = (max(max(simMx))==1) | sum(basal)~=nBasal;
        end
    
        
        B0 = .95*rand(S,1)+.05;
        
        SList = 1:S;
        idxPar = datasample(SList(~basal),sum(~basal),'Replace',false);
        
        %short-weighted trophic level
        A = full(sparse(res,con,1,S,S));
        nonCannA = A;
        nonCannA(diag(true(S,1)))=0;

        %Formula for the prey-averaged trophic level.
        patl_mx = sparse(nonCannA)*(diag(1./sum(sparse(nonCannA))));
        patl = (speye(S)-patl_mx')\ones(S,1);
        
        webData{ii}.web = ii;
        webData{ii}.B0 = B0;
        webData{ii}.LL = [res,con];
        webData{ii}.gr = gr;
        webData{ii}.patl = patl;
        webData{ii}.parOrder = idxPar;
    end
    
    for  model = allModels'
        
        kFree = kFrees(model(1));
        kPara = kParas(model(2));
        
        fracFree = fracFrees(model(3));
        fracPara = fracParas(model(4));

        for fPar = fParAll0
            %We don't do simulations we have already done
            if sum(kFreesDone==kFree)&&sum(kParasDone==kPara)&&sum(fParDone==fPar)&&(ii<=nWebsDone)
                continue
            end
            %Most factors don't matter if fPar==0.
            if (fPar==0)&&((model(3)~=1)||(model(4)~=1)||(model(2)~=1))
                continue
            end
            
            %Number of parasites at this level. went with round to keep as
            %close to what was advertised as possible.
            nBasal = sum(basal);
            nFree = S-nBasal;
            para = false(S,1);
            nPar = round(fPar*nFree);
            para(idxPar(1:nPar)) = true;
            %simParams is what gets split up by the parfor loop and defines how each simulation goes.
            simNo = simNo+1;
            simParams{simNo}.web = ii;
            simParams{simNo}.fPar = fPar;
            simParams{simNo}.kFree = kFree;
            simParams{simNo}.kPara = kPara;
            simParams{simNo}.fracFree = fracFree;
            simParams{simNo}.fracPara = fracPara;
            simParams{simNo}.modelCode = [fracFree fracPara];
            simParams{simNo}.para = para;
            simParams{simNo}.LL = [res con];
            simParams{simNo}.B0 = B0;
            simParams{simNo}.gr = gr;
            simParams{simNo}.parOrder = idxPar;
            simParams{simNo}.patl = patl;
            
        end
        
    end
end

%Save the webData.
save('../raw/webData.mat','webData')

%%%All these parameters are the same for all simulations. Can change these, but it would have to be done
%%%manually (i.e. would need to copy webData files to a new directory, modify code so that it loads
%%%those webs and runs all simulations on those webs with the new parameters... yeesh.)

%Metabolic parameters for body sizes, metabolic rates
axFree = .314;
axPara = .314;
metabScaling = 0.25;

%%% Parameters needed in dynamical equations
%maximum production rate of consumers relative to metabolic rate
yFree = 8;
yPara = 8;

%shared carrying capacity for basal species
K = 5;

%Hill coefficient and half-saturation density for the functional response
h=1.2;
halfSat = 0.5;

%Simulation time is: Max(Tf,t_(last Extinction)+Trelax)
Tf = 10000;
Trelax = 2000;

%%% Simulation parameters
%Specifying threshold for species extinction and tolerance for solvers.
%These are necessarily related; must have abstol smaller that extctThresh
%so that we can have confidence in the species actually going extinct;
%otherwise we are guaranteed no accuracy in the species that goes extinct!
%Sadly this can have a pretty significant impact on run time.
extctThresh = 1e-10;
AbsTol = 1e-16;
RelTol = 1e-6;

%This is the structure that passes into the integrate parasites function.
%contains all parameters needed to run the simulation.

params = struct(... 
    ... Web Parameters:
    'S',S...
    ,'C',C...
    ... DS Parameters:
    ,'K',K...
    ...FR Parameters:
    ,'halfSat',halfSat...
    ,'phi',.15...
    ,'h',h...
    ,'eijBas',0.45 ...
    ,'eijCon',0.85 ...
    ,'yFree',8 ...
    ,'yPara',8 ...
    ,'axFree',axFree ...
    ,'axPara',axPara ...
    ... Link Level properties:
    ,'res',[]...
    ,'con',[]...
    ,'eij',[]...
    ,'wij',[] ...
    ,'yij',[]...
    ... Species level Properties:
    ,'basal',zeros(40,1)...
    ,'r',zeros(40,1)...
    ,'B0',zeros(40,1)...
    ,'x',zeros(40,1)...
    ,'M',zeros(40,1)...
    ,'modelCode',[0 0]...
    ,'para',zeros(40,1)...
    ... Simulation Parameter:
    ,'Tf',Tf ...
    ,'Trelax',Trelax ...
    ... Solver Parameters:
    ,'extctThresh',extctThresh ...
    ,'AbsTol',AbsTol...
    ,'RelTol',RelTol...
    ,'odeSolver',@ode45... .
    ,'options',odeset()...
    );

[TS, extcts,bs,rs] = integrateExperiment(simParams,params,S,nSims);

if appending %load old outputs and add the new ones.
    oldSimParams = load('../raw/simParams.mat');
    simParams = [oldSimParams;simParams];

    oldData = load('../raw/rawOutputs.mat');
    TS = cat(3,oldData.TS,TS);
    extcts = [oldData.extcts; extcts];
    bs = [oldData.bs;bs];
    rs = [oldData.rs;rs];
end

%Save the outputs
save('../raw/simParams.mat','simParams');
save('../raw/rawOutputs.mat','extcts','TS','bs','rs','-v7.3');

%Save quite a bit of time by running these abundance curves in parallel.

%process the data to a form that is a bit more manageable.
run('parseData.m')

%save the plotting data.
run('makePlottingData.m')

%save the new data.
run('fullBreakdown.m')
