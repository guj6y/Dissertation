%It isn't at all clear what you load here. So:
%out.mat contains extct and TS.  %  extct is a cell array with nSims cells. Each cell contains a 41xn array. (1,:) slices the event times (including start and end). (2:end,ii) slices the biomasses of all species at the ii-th extinction event.
%  TS is a 3-d array of size 40x1000xnSims. This is the last 1000 time steps of the biomasses for all species. This should be useful for analyzing the equilibrium state (check for oscillations, fft,etc.)
%
%simParams
%.mat contains the cell array simParams. each element of the cell array is a structure. This is the structure:
%{
          web: 1                 %This is the index of tyhe niche web used
         fPar: 0                 %This is the fraction of parasites 
        kFree: 1                 %This is exponent  on free-living BSR
        kPara: -3                %This is exponent on parasite BSR
     fracFree: 0                 %This is host refuge identifier   
     fracPara: 0                 %This is the concomittant identifier
    modelCode: [0 0]             %This is redundant.
         para: [40×1 logical]    %This tells which species are parasites
           LL: [250×2 double]    %This is the link list
           B0: [40×1 double]     %THis is the initial Biomass.
           gr: [40×1 double]     %Growth rates of basal species
     parOrder: [1×33 double]     %This is the order of parasites.
           TL: [40×1 double]     %This is nothing.
         patl: [40×1 double]     %This is the prey-averaged trophic level.
%}
if ~exist('simParams','var')
    load '../raw/rawOutputs.mat'
    load '../raw/simParams.mat'
    load '../raw/metaSimData.mat'
end

%The goal of this code is to convert the contents of rawoutputs.mat to complicated 
%multi-dimensional arrays that I can plot. This code should generate all the data
%that will eventually be plotted.
nFPar = numel(fParAll0);
nSims = numel(extcts);

nanArray = nan(nWeb,nFPar,nFacts(1),nFacts(2),nFacts(3),nFacts(4));

abc = struct('common',struct('slopes',nanArray...
                            ,'RSquared',nanArray)...
            ,'natural',struct('slopes',nanArray...
                             ,'rSquared',nanArray)...
            );


nanStruct = struct('all',nanArray...
                  ,'free',nanArray...
                  ,'para',nanArray...
                  ,'basal',nanArray...
                  );

persistences = nanStruct;
cvs = nanStruct;
biomasses = nanStruct;

frac = struct('free',nanArray,...
                            'para',nanArray,...
                            'basal',nanArray);
                             
activities = nanStruct;
activities.con = nanArray;

nanSpeciesArray = nan(S,nWeb,nFPar,nFacts(1),nFacts(2),nFacts(3),nFacts(4));

speciesTypes = struct('para',nanSpeciesArray...
                     ,'free',nanSpeciesArray...
                     ,'basal',nanSpeciesArray...
                     );

simNumbers = nanArray;

for ii = 1:nSims
    webNo = simParams{ii}.web;
    fact1Level = simParams{ii}.kFree == kFrees;
    fact2Level = simParams{ii}.kPara == kParas;
    fact3Level = simParams{ii}.fracFree == fracFrees;
    fact4Level = simParams{ii}.fracPara == fracParas;
    fParLevel = simParams{ii}.fPar == fParAll0;
    
    thisSim_web = {webNo,fParLevel,fact1Level,fact2Level,fact3Level,fact4Level};
    thisSim_S = [1:S,thisSim_web];

    simNumbers(thisSim_web{:}) = ii; 
    para = simParams{ii}.para;
    basal = simParams{ii}.gr>0;
    free = ~(para|basal);
    
    M = zeros(S,1);
    x = zeros(S,1);
    kFree = kFrees(fact1Level);
    kPara = kParas(fact2Level);
    patl = simParams{ii}.patl;
    M(free) = (10.^kFree).^(patl(free)-1);
    M(para) = 10.^(kPara + kFree*(patl(para)-2));
    x = .314.*M.^(-0.25);
    x(basal) = simParams{ii}.gr(basal);
     
    speciesTypes.para(thisSim_S{:})  = para;
    speciesTypes.free(thisSim_S{:})  = free;
    speciesTypes.basal(thisSim_S{:}) = basal;
    
    abc.common.slopes(thisSim_web{:})  = bs(ii,1);
    abc.natural.slopes(thisSim_web{:}) = bs(ii,2);

    abc.common.rSquared(thisSim_web{:})  = rs(ii,1);
    abc.natural.rSquared(thisSim_web{:}) = rs(ii,2);


    
    %Total biomass is weird to me. Why? -glad you asked. It's because the scale of species biomasses
    %varies geometrically. So simply summing up to get a total biomass ends up being skewed by a 
    %few very abundant species, with the last few decimals changing slightly because of most of 
    %the rest - the rest that have a profound impact on the overall state we are in. AT least, they
    %have a profound impact on the persistence metric! 
    meanBiomasses = mean(TS(:,:,ii),2);
   
    persistences.all(thisSim_web{:})   = mean(meanBiomasses>0);
    persistences.para(thisSim_web{:})  = mean(meanBiomasses(para)>0);
    persistences.free(thisSim_web{:})  = mean(meanBiomasses(free)>0);
    persistences.basal(thisSim_web{:}) = mean(meanBiomasses(basal)>0);

    numParasites   = sum(meanBiomasses(para) > 0);
    numFreeLivers = sum(meanBiomasses(free) > 0);
    numBasal = sum(meanBiomasses(basal) > 0);
    numAlive = sum(meanBiomasses > 0);
    
    frac.para(thisSim_web{:}) = numParasites./numAlive;
    frac.free(thisSim_web{:}) = numFreeLivers./numAlive;
    frac.basal(thisSim_web{:}) = numBasal./numAlive;

    allBiomassTS   = sum(TS(:,:,ii),1);
    paraBiomassTS  = sum(TS(para,:,ii),1);
    freeBiomassTS  = sum(TS(free,:,ii),1);
    basalBiomassTS = sum(TS(basal,:,ii),1);
    
    meanAllBiomass   = mean(allBiomassTS);
    meanParaBiomass  = mean(paraBiomassTS);
    meanFreeBiomass  = mean(freeBiomassTS);
    meanBasalBiomass = mean(basalBiomassTS);
    
    stdAllBiomass   = std(allBiomassTS);
    stdParaBiomass  = std(paraBiomassTS);
    stdFreeBiomass  = std(freeBiomassTS);
    stdBasalBiomass = std(basalBiomassTS);

    cvAllBiomass   = stdAllBiomass./meanAllBiomass;
    cvParaBiomass  = stdParaBiomass./meanParaBiomass;
    cvFreeBiomass  = stdFreeBiomass./meanFreeBiomass;
    cvBasalBiomass = stdBasalBiomass./meanBasalBiomass;

    biomasses.all(thisSim_web{:})   = meanAllBiomass;
    biomasses.para(thisSim_web{:})  = meanParaBiomass;
    biomasses.free(thisSim_web{:})  = meanFreeBiomass;
    biomasses.basal(thisSim_web{:}) = meanBasalBiomass;

    cvs.all(thisSim_web{:})   = cvAllBiomass;
    cvs.para(thisSim_web{:})  = cvParaBiomass;
    cvs.free(thisSim_web{:})  = cvFreeBiomass;
    cvs.basal(thisSim_web{:}) = cvBasalBiomass;
    
    freeActivityTS  = sum(TS(free,:,ii).*x(free),1);
    paraActivityTS  = sum(TS(para,:,ii).*x(para),1);
    basalActivityTS = sum(TS(basal,:,ii).*x(basal),1).*(1-basalBiomassTS/5);
    conActivityTS   = freeActivityTS + paraActivityTS;
    totalActivityTS = conActivityTS + basalActivityTS;
    
    meanAllActivity   = mean(totalActivityTS);
    meanParaActivity  = mean(paraActivityTS);
    meanFreeActivity  = mean(freeActivityTS);
    meanBasalActivity = mean(basalActivityTS);
    meanConActivity   = mean(conActivityTS);
    
    activities.all(thisSim_web{:})   = meanAllActivity;
    activities.para(thisSim_web{:})  = meanParaActivity;
    activities.free(thisSim_web{:})  = meanFreeActivity;
    activities.basal(thisSim_web{:}) = meanBasalActivity;
    activities.con(thisSim_web{:})   = meanConActivity;

end

save('../raw/out.mat','persistences','biomasses','cvs','speciesTypes','abc','activities');
