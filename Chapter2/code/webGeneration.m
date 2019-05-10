#This script transforms data from http://esapubs.org/archive/ecol/... and
#creates food webs (networks of trophic interactions between species).  I
#want to consistenly aggregate and select the various species.  It saves
#the created webs as link-lists:
#
# res,con,para
#  i,j,b
#
#Here, i is the resource index, j is the consumer index, and b is a binary
#variable that is 1 if the link is parasitic (i not a parasite, j is a
#parasite).  I also save the web properties from calculateLocalProperites.

#I had to do quite a bit of post-processing on the files inorder to make
#them somewhat consistent and not a nightmare to import (the
#post-processing was a nightmare in itself, but whatevs)


#Detailed descriptions of the column headings and variables are available
#in the metadata.

#Requires a somewhat recent
clear
linkFiles = {'BSQweb_Links.txt',...
    'CSMweb_Links.txt',...
    'EPBweb_Links.txt',...
    'Flensburg_Data_Links.txt',...
    'Otago_Data_Links.txt',...
    'Sylt_Data_Links.txt',...
    'Metaweb_Links.csv'};

numColLinks = [20,20,20,16,16,16,12];

nodeFiles = {'BSQweb_Nodes.txt',... Has body sizes
    'CSMweb_Nodes.txt',...          Has body sizes
    'EPBweb_Nodes.txt',...          Has body sizes
    'Flensburg_Data_Nodes.txt',...  No body sizes
    'Otago_Data_Nodes.txt',...      NO body sizes
    'Sylt_Data_Nodes.txt',...
    'Metaweb_Nodes.csv'};        #No body sizes

numColNodes = [46,46,46,41,41,41,46];

OGConBSCell = cell(3,1);
OGResBSCell = cell(3,1);
OGLinkTypeCell = cell(3,1);

nWebs = 6;

if ~exist('concomittantWebs','var')
    concomittantWebs = false;
end

if concomittantWebs 
    #This is probably better. Those concomittant links don't belong in
    #a web of feeding relationships.
    fprintf('Generating webs with concomittant links\n');
else
    fprintf('Generating webs without concomittant links\n');
end
topoConcomittant = false
codeDir = pwd;

shortWebNames = {'bahia','carp','punta','flens','otago','sylt','meta'};

#These are the links that I include in the food web.  We skip non-trophic
#interactions, but it might be important to double-check which links I use.
# Descriptions are found in the metadata.
linkTypesTrophic = [1 ...   predation
    ,2 ...                  Social Predation
    ,3 ...                  MicroPredation
    ,18 ...                 Detritivory
    ,21 ...                 Facultative Micropredation
    ];

linkTypesTrophicPara = [16 ...  Predation on free-living non-feeding stage
    ,19 ...                     Parasite Intraguild Antagonism
    ];

linkTypesPara = [4 ...  Parasitic Castration
    ,5 ...              Pathogen Infection
    ,6 ...              Macroparasitism
    ,8 ...              Pollination
    ,10 ...             Trophically Transmitted Parasitic Castration
    ,11 ...             Trophically Transmitted pathogen infection
    ,12 ...             Trophically Transmitted Parasitism
    ];

concLinks = [14 ...     Concurrent predation on symbionts
    ,15 ...             Trophic Transmission
    ];

if concomittantWebs
    linkTypesIncluded = [linkTypesTrophic linkTypesTrophicPara...
        linkTypesPara concLinks];
else
    linkTypesIncluded = [linkTypesTrophic linkTypesTrophicPara...
        linkTypesPara];
end



propertiesCell = cell(nWebs,2);
oldNamesCell = cell(nWebs,1);
newNamesCell = cell(nWebs,1);

oldClassCell = cell(nWebs,1);
newClassCell = cell(nWebs,1);

nEndo = zeros(nWebs,1);
totSpecies = zeros(nWebs,1);
totPara = zeros(nWebs,1);
totBasal = zeros(nWebs,1);
totFreeCon = zeros(nWebs,1);

Ls = zeros(nWebs,1);
connectances = zeros(nWebs,1);
binEndoCell = cell(nWebs,1);
binEctoInvertCell = cell(nWebs,1);
binEctoVertCell = cell(nWebs,1);
binAutoCell = cell(nWebs,1);
binParaCell = cell(nWebs,1);

binOldEndoCell = cell(nWebs,1);
binOldEctoInvertCell = cell(nWebs,1);
binOldEctoVertCell = cell(nWebs,1);
binOldAutoCell = cell(nWebs,1);
binOldDetCell = cell(nWebs,1);
binOldParaCell = cell(nWebs,1);

speciesTypeCell = cell(nWebs,1);
meanEndoPara = zeros(nWebs,1);

nodeDatas = cell(nWebs,1);
linkDatas = cell(nWebs,1);
linkListCell = cell(nWebs,1);
aggregationLevel = 1;
jjLinksNInitial = zeros(nWebs,length(linkTypesIncluded)+1);
linkTypesCell = cell(nWebs,1);
linkParaCell = cell(nWebs,1);

meanGen = zeros(nWebs,2);
varGen = zeros(nWebs,2);
seGenPool = zeros(nWebs,1);

OGSpeciesLinks = cell(nWebs,1);

for ii = 1:nWebs
    cd('data/Raw/')
    cd('Links')
    if ii<7
        # linkData = readtable(linkFiles{ii},'Delimiter','tab','ReadVariableNames',true);
        # No readtable in octave; Relatively quick workaround that preserves the
        # rest of the code could be to create a structure with column headers.
        link_fid = fopen(linkFiles{ii});
        link_headers = textscan(link_fid, '%s', 1, 'Delimiter', '');
        link_headers = link_headers{1}{1};
        link_headers = textscan(link_headers,'%s',-1,'Delimiter','\t');
        link_headers = link_headers{1};
        link_headers = cellfun(@(x) strrep(x,' ',''), link_headers, 'UniformOutput',false);
        link_data = textscan(link_fid,'%s','Delimiter','');
        link_data = cellfun(@(x) textscan(x,'%s',...
                                          'Delimiter','\t',...
                                          'MultipleDelimsAsOne',false)...
                                          ,link_data{1});
        link_data = [link_data{:}]';
        
        
    else
        # linkData = readtable(linkFiles{ii},'Delimiter','comma','ReadVariableNames',true);
    end
    if (ii>3)&&(ii<7)
        #linkData.Properties.VariableNames(7) = {'LinkTypeID'};
        link_headers{7} = 'LinkTypeID';
        
    end
    linkData = link_data;
    linkDatas{ii} = linkData;
    #linksIncluded is a binary array that indicates the links to include
    #in the initial food web.
    col_dict_link = containers.Map(link_headers,1:numel(link_headers));
    linksIncluded = zeros(size(linkData(:,col_dict_link('LinkTypeID'))));

    #Detritus?

    #These Nodes are detritus - I need to be sure to take them out
    #completely since I was getting detritus as a basal species (!?) - some
    #links of snails eating detritus were labeled as predation; 4 webs have
    #detritus, 3 label 'detritus', the  other as 'Detritus'
    #detritusNodes = strcmp(nodeData.OrganismalGroup,'Detritus')|...
    #    strcmp(nodeData.OrganismalGroup,'detritus');

    # detritusNodeIDs = nodeData.NodeID(detritusNodes);
    # detritusSpeciesIDs = speciesID(detritusNodes);

    #linksPara is a binary array that indicates whether the link type
    #implies a parasitic consuemr
    linksPara = linksIncluded;
    linksConc = linksIncluded;
    #Extract the links that I want; see metadata.
    count = 0;
    link_type_array = cellfun(@(x) str2num(x),linkData(:,col_dict_link('LinkTypeID')));
    for jj = linkTypesIncluded
        count = count+1;
        jjLinks = link_type_array == jj;
        linksIncluded = (linksIncluded) + (jjLinks)*jj;
        jjLinksNInitial(ii,count) = sum(jjLinks);
        if sum(linkTypesPara==jj)
            linksPara = (linksPara) | (jjLinks);
        end
        if jj == 14
            linksConc = linksConc | jjLinks;
        end

    end
    jjLinksNInitial(ii,end) = sum(jjLinksNInitial(ii,1:(end-1)),2);
    #List of IDs of resources, consumers, and parasites as species; res0
    #and con0 define the foodweb; para0 will become a parasitic species
    #identifier.
    res0 = linkData(:,col_dict_link('ResourceSpeciesID'))(linksIncluded>0);
    con0 = linkData(:,col_dict_link('ConsumerSpeciesID'))(linksIncluded>0);
    para0 = linkData(:,col_dict_link('ConsumerSpeciesID'))(linksPara);

    OGSpeciesLinks{ii} = {[res0,con0],para0};
    #Trying to add all conceivable concomittant links.  This tries to add
    #about 16K links for the first web... which is a *bit* high.  Don't do
    #this?
    if topoConcomittant
        resToAdd = [];
        conToAdd = [];
        for jj = unique(para0')
            jjsHosts = res0(con0==jj);
            if numel(jjsHosts >0)
                for kk = jjsHosts'
                    kksPred = con0(res0==kk);
                    concAdded = [];
                    for ll = kksPred'
                        if sum(ll==jj)>0
                            continue
                        elseif sum(ll==jjsHosts) >0
                            continue
                        else
                            concAdded = [concAdded; ll];
                        end
                    end
                    resToAdd = [resToAdd; ones(size(concAdded))*jj];
                    conToAdd = [conToAdd; concAdded];
                end
            end

        end
        res0 = [res0; resToAdd];
        con0 = [con0; conToAdd];
    end
    res0 = cellfun(@(x) str2num(x), res0);
    con0 = cellfun(@(x) str2num(x), con0);
    #Identifying Concomittant links; still not sure about intermediate host
    #at this point.  Before the following two lines, these are the same
    #size as the original link list from data.
    linksConc = linksConc(linksIncluded>0);
    linksPara = linksPara(linksIncluded>0);

    #Now, linksConc is a binary array of the same length as the initial
    #link list (i.e. it is L x 1 vector), =1 if the link is concomittant,
    #and zero otherwise.  Similarly for linksPara.

    #linkIDs Specifies the type of each link.  Could very well use the
    #original link types from the data, too!  Do this for now since I'm not
    #overly concerned about all those different types.
    linkTypes = linksIncluded(linksIncluded>0);
    linkIDs =  ones(size(linksConc)) + linksPara + 2*linksConc;

    #Actually Removing detrital links
    #for jj = detritusSpeciesIDs'
    #
    #    jjDetritalLinks = (res0==jj)|(con0==jj);
    #    res0(jjDetritalLinks) = [];
    #    con0(jjDetritalLinks) = [];
    #    para0(para0==jj) = [];

    #end


    #Unique IDs of resource, consumers, as stages:
    #resStage0 = linkData.ResourceNodeID(linksIncluded);
    #conStage0 = linkData.ConsumerNodeID(linksIncluded);
    #paraStage0 = linkData.ConsumerNodeID(linksPara);
    #
    #Actually Removing Detrital links ( and therefore species)
    #for jj = detritusNodeIDs'
    #    jjDetritalLinks = (res0==jj)|(con0==jj);
    #    resStage0(jjDetritalLinks) = [];
    #    conStage0(jjDetritalLinks) = [];
    #    paraStage0(paraStage0==jj) = [];
    #
    #end

    #Read node data
    
    cd('../Nodes')
    if ii<7
        node_fid = fopen(nodeFiles{ii});
        node_headers = textscan(node_fid, '%s', 1, 'Delimiter', '');
        node_headers = node_headers{1}{1};
        node_headers = textscan(node_headers,'%s',-1,'Delimiter','\t');
        node_headers = node_headers{1};
        node_headers = cellfun(@(x) strrep(x,' ',''), node_headers, 'UniformOutput',false);
        node_data = textscan(node_fid,'%s','Delimiter','');
        node_data = cellfun(@(x) textscan(x,'%s',...
                                          'Delimiter','\t',...
                                          'MultipleDelimsAsOne',false)...
                                          ,node_data{1});
        node_data = [node_data{:}]';
        # nodeData = readtable(nodeFiles{ii},'Delimiter','tab','ReadVariableNames',true);
    else
        nodeData = readtable(nodeFiles{ii},'Delimiter','comma','ReadVariableNames',true);
    end
    cd('..');
    col_dict_node = containers.Map(node_headers,1:numel(node_headers));
    nodeData = node_data;
    nodeDatas{ii} = nodeData;
    node_ids = nodeData(:,col_dict_node('NodeID'));
    [~, I] = sort(node_ids);
    nodeData = nodeData(I,:);
    
    
    
    # I hate using unique, but..
    idList = unique([res0;con0]);
    basalIds = setdiff(res0,con0);
    
    bodySizes = zeros(1,max(idList));
    biomasses = zeros(1,max(idList));
    
    speciesID = nodeData(:, col_dict_node('SpeciesID'));
    speciesID = cellfun(@(x) str2num(x), speciesID);
    stageID = nodeData(:, col_dict_node('StageID'));
    stageID = cellfun(@(x) str2num(x), stageID);
   
    # Doing the bodysizes here.
    if ii < 4
      bodySizes_g = nodeData(:, col_dict_node('BodySize(g)'));
      bodySizes_g = cellfun(@(x) str2num(strjoin({'0',x},'')), bodySizes_g);
      biomass_kg_ha = nodeData(:, col_dict_node('Biomass(kg/ha)'));
      biomass_kg_ha = cellfun(@(x) str2num(strjoin({'0',x},'')), biomass_kg_ha);
      for jj = idList'
          #The body size of species j is the adult body size of species j.
          #The biomass of species j is the sum of the biomasses of all
          #stages.  NaN are omitted; biomass not that important, mainly as a
          #decent initial condition for simulations on these webs. If the
          #adult does not have a body size, take the body size of the most
          #abundant stage. If the most abundant stage does not have a body
          #size, take the lowest stage IDs body size.  Hopefully the only
          #thing that this leaves as NaN are the few messed up species, and
          #the plants.

          binj = speciesID == jj;
          minStage = min(stageID(binj));
          binMin = stageID == minStage;

          bodySizes(jj) = mean(bodySizes_g((binj)&(binMin)));

          biomasses(jj) = nansum(biomass_kg_ha(binj));
          if sum(binj) > 1
              if isnan(bodySizes(jj))&&sum(isfinite(biomass_kg_ha(binj))>0)&&~(sum(isnan(nodeData.BodySize_g_((binj)&(binMin))))==sum((binj)&(binMin)))
                  bodySizes(jj) = bodySize_g((binj)&...
                      (biomass_kg_ha==...
                      max(biomass_kg_ha(binj&(isfinite(bodySize_g))))...
                      )&isfinite(bodySize_g));
              end


          end
      endfor
    endif
    
    #Some years later, came back andthought it woud be nice to look at
    #OG data. this turns out to be a much simpler way.
    node_class = nodeData(:, col_dict_node('Class'));
    if ii<=3
    nNodes = numel(node_ids);
    nodeToIdx = containers.Map(nodeData(:, col_dict_node('NodeID')),1:nNodes);
    
    OGCon = linkData(:,col_dict_link('ConsumerNodeID'));
    OGRes = linkData(:,col_dict_link('ResourceNodeID'));
    OGLinkType = linkData(:,col_dict_link('LinkTypeID'));
    OGLinkType = cellfun(@(x) str2num(x), OGLinkType);
    OGGoodLinks = sum(OGLinkType == linkTypesIncluded,2)>0;
    OGParaLinks = sum(OGLinkType == linkTypesPara,2)>0;
    
    OGConBS = bodySizes_g(cellfun(@(x) nodeToIdx(x),OGCon));
    OGResBS = bodySizes_g(cellfun(@(x) nodeToIdx(x),OGRes));
    
    OGEctoVert = sum([cellfun(@(x)~isempty(x),strfind(node_class,'Agnatha'))...
                     ,cellfun(@(x)~isempty(x),strfind(node_class,'Chondrichthyes'))...
                     ,cellfun(@(x)~isempty(x),strfind(node_class,'Ostheichthyes'))...
                    ,cellfun(@(x)~isempty(x),strfind(node_class,'Amphibia'))...
                    ,cellfun(@(x)~isempty(x),strfind(node_class,'Reptilia'))...
                    ,cellfun(@(x)~isempty(x),strfind(node_class,'Osteichthyes'))...
                  ],2);
    
    OGEndo = sum([cellfun(@(x)~isempty(x),strfind(node_class,'Aves'))...
                     ,cellfun(@(x)~isempty(x),strfind(node_class,'Mammalia'))...
                     ],2);
    
    OGPara = false(nNodes,1);
    OGPara( cellfun(@(x) nodeToIdx(x),OGCon(OGParaLinks))) = true;
                 
    OGInvert = ~(OGEctoVert|OGEndo|OGPara);
    
    OGLinkType = OGInvert + 2*OGEctoVert + 3*OGEndo + 4*OGPara;
    OGLinkType = OGLinkType(cellfun(@(x) nodeToIdx(x),OGCon));
    
    OGConBS = OGConBS(OGGoodLinks);
    OGResBS = OGResBS(OGGoodLinks);
    OGLinkType = OGLinkType(OGGoodLinks);
    definedBS = isfinite(OGConBS)&isfinite(OGResBS);
    OGConBSCell{ii} = OGConBS(definedBS);
    OGResBSCell{ii} = OGResBS(definedBS);
    OGLinkTypeCell{ii} = OGLinkType(definedBS);
    end
    #The total number of species in the Data set; the final number included
    #in my foodweb may be less than this!!  Actually it's not excactly
    #that for carpinteria, punta and bahia because of how their data is
    #saved
    nAllSpecies = max(speciesID);

    #SEtting up names
    oldNamesCell{ii} = cell(nAllSpecies,1);
    working_name = nodeData(:, col_dict_node('WorkingName'));
    oldNamesCell{ii}(speciesID) = working_name;

    oldClassCell{ii} = cell(nAllSpecies,1);
    oldClassCell{ii}(speciesID) = node_class;

    
    binOldEndoCell{ii} = zeros(nAllSpecies,1);
    binOldEctoVertCell{ii} = zeros(nAllSpecies,1);
    binOldEctoInvertCell{ii} = zeros(nAllSpecies,1);
    binOldAutoCell{ii} = zeros(nAllSpecies,1);
    binOldDetCell{ii} = zeros(nAllSpecies,1);

    [oldClassCell{ii}{~cellfun(@ischar,oldClassCell{ii})}] = deal('');
    #Ectotherm Vertebrates: one of these classes
    binAgna = ~cellfun('isempty',strfind(oldClassCell{ii},'Agnatha'));
    binChon = ~cellfun('isempty',strfind(oldClassCell{ii},'Chondrichthyes'));
    binOsth = ~cellfun('isempty',strfind(oldClassCell{ii},'Ostheichthyes'));
    binOste = ~cellfun('isempty',strfind(oldClassCell{ii},'Osteichthyes'));
    binAmph = ~cellfun('isempty',strfind(oldClassCell{ii},'Amphibia'));
    binRept = ~cellfun('isempty',strfind(oldClassCell{ii},'Reptilia'));
    binOldEctoVertCell{ii} = binAgna|binChon|binOste|binAmph|binRept|binOsth;

    #Endotherms: one of these classes
    binAvesCell = ~cellfun('isempty',strfind(oldClassCell{ii},'Aves'));
    binMammCell = ~cellfun('isempty',strfind(oldClassCell{ii},'Mammalia'));
    binOldEndoCell{ii} = binAvesCell|binMammCell;
    #sum(binEndoCell{ii}&binEctoVertCell{ii})
    #Autotrophs: generality = 0
    node_feeding = nodeData(:, col_dict_node('Feeding'));
    node_org_group = nodeData(:, col_dict_node('OrganismalGroup'));
    binOldAutoCell{ii}(speciesID(~cellfun('isempty',strfind(node_feeding,'autotroph'))))=1;
    binOldDetCell{ii}(speciesID(~cellfun('isempty',strfind(node_org_group,'detritus'))))=1;

    #sum(binEndoCell{ii}&binAutoCell{ii})
    #sum(binEctoVertCell{ii}&binAutoCell{ii})
    #Ectotherm invertebrates: the rest.  There are tiny little critters
    #that maybe shouldn't be classified as ectotherm invertebrates (i.e.
    #single cell organisms).. there aren't very many nodes, but might they
    #be important enough to delineate?
    binOldEctoInvertCell{ii} = ...
        ~(binOldEctoVertCell{ii}|binOldEndoCell{ii}|binOldAutoCell{ii}|binOldDetCell{ii});
    
    #oldNamesCell{ii}(~binSpeciesIncluded) = [];
    #List of all species Ids in the datase
    allSpeciesList = 1:nAllSpecies;
    
    #Binary array that indicates whether a particular species out of the
    #entire data set is included in my food web.  #WOrking with species Ids
    binSpeciesIncluded = false(nAllSpecies,1);
    
    ectoOnly=false;
    if ectoOnly
        noEcto = ~((binOldEctoInvertCell{ii}(res0)&binOldEctoInvertCell{ii}(con0))...
            |(binOldAutoCell{ii}(res0)&binOldEctoInvertCell{ii}(con0))...
            |(binOldDetCell{ii}(res0)&binOldEctoInvertCell{ii}(con0))...
            );
        res0(noEcto)=[];
        con0(noEcto)=[];
    end
    binSpeciesIncluded([res0;con0]) = true;
    biomasses = biomasses(binSpeciesIncluded);
    bodySizes = bodySizes(binSpeciesIncluded);
    #List of only the old species Ids that are included in my food web.
    oldIdsInWeb = allSpeciesList(binSpeciesIncluded);
    #Number of species in new web:
    nSpecies = length(oldIdsInWeb);
    #List of new IDs in web
    newIds = 1:nSpecies;

    #New names in web
    newNamesCell{ii} = cell(nSpecies,1);
    newClassCell{ii} = cell(nSpecies,1);

    #The replacement vector works by placing the new ID at the old ID's
    #position.
    replacementVector = zeros(nAllSpecies,1);
    replacementVector(oldIdsInWeb) = newIds;

    #replacementVector('oldID of species j') = 'new ID of species j'
    #so,this re-does all the IDs correctly.  This doesn't change the length
    #of res or con, so now the link list has redundancies.  It doesn't
    #really matter that much since we can use 'sparse' to get the unique
    #link list.
    res = replacementVector(res0);
    con = replacementVector(con0);
    
    noLink = (res==0)|(con==0);
    res=res(~noLink);
    con=con(~noLink);
    #Note that we don't need to worry link IDs since the lengths aren't
    #changing, so we still have the links correctly identified.  However,
    #there is still the question of whether or not concomittant links get
    #grouped with non-concomittant links.


    #Doesn't work so nicely for names,sizes, and biomasses.
    newNamesCell{ii} = oldNamesCell{ii}(oldIdsInWeb);
    newClassCell{ii} = oldClassCell{ii}(oldIdsInWeb);

    #Binary array that indicates whether a particular species is a
    #parasite; out of all species in the dataset
    binParaAll = false(nAllSpecies,1);
    para0 = cellfun(@(x) str2num(x), para0);
    binParaAll(para0) = true;

    #From that, only taking the species that are actually included.
    para = binParaAll(binSpeciesIncluded);

    #Is the final web a *trophic* web?
    trophicWeb = false;
    [LL, idxLL] = sortrows([con res]);
    con = LL(:,1);
    res = LL(:,2);
    linkTypes = linkTypes(idxLL);
    cd(codeDir)
    #Until the web is trophic
    while ~trophicWeb
        #CAUTION: this could combine parasites with free-livers.
        #Probably not something we want to do.  It doesn't happen in the
        #six original ecosystems I looked at (punta, bahia, carpinteria,
        #sylt, flensburg, otago), but still a concern if this code gets
        #recycled.  I highly doubt it would happen, but you never know...
        
        #TODO: Apropos of above, stop this code from combining 
        #parasites with free-livers. 
        
        #Is that ^ not an indication that they are a functionally distinct
        #class of species?
        
        #FInd disconnected species: This should only trigger once...
        S_ = max([res;con]);
        A = sparse(res,con,1,S_,S_);
        # [n,grps] = graphconncomp(sparse(A),'weak',true,'directed',true);
        [GC, I] = giantComponent(sparse(A+A'));
        #{
        maxCompSize = 0;
        compSizes = zeros(1,n);
        for ll = 1:n
            sizeii = sum(grps==ll);
            compSizes(ll) = sizeii;
            if sizeii>maxCompSize
                maxCompSize=sizeii;
                largeComp = ll;
            end
        end
        #}
        winners = false(S_,1);
        winners(I) = true;
        [killTheseI, killTheseJ] = find((1-winners)'*(1-winners));
        #similarity Matrix
        simMx = calculateSimilarity(res,con);
        #Same matrix:  (sim(i,j) = 1, if i is the same as j.& =0 o.w.
        sameMx = triu(simMx>=aggregationLevel,1);
        
        #we want I -> J; I lists indices of superfluous species.
        [I,J] = find(sameMx);
        [S1,S2] = size(sameMx);
        badCombo = (para(I) == ~para(J)) | (~para(I) == para(J));
        I(badCombo) = [];
        J(badCombo) = [];
        
        sameMx = sparse(I,J,1,S1,S2);
        
        if (sum(sum(sameMx)))==0
            trophicWeb = true;
            continue
        end
        #number of species in the web.
        nSpecies = max([res;con]);

        #Species that we need to delete; the triu means that doing this
        #won't delete any species equivalence classes (trophic species)
        deleteVector = (1:nSpecies)';
        I = [I;killTheseI];
        J = [J;killTheseJ];
        deleteVector(I) = 0;

        #Delete links that include a superfluous species; first, change
        #the species to zero; res0 and con0 are 1 if the link is to be
        #deleted, and zero otherweise.
        res0 = deleteVector(res)==0;
        con0 = deleteVector(con)==0;


        #Then, delete all links that contain a zero; must do this to both
        #the resource and consumer vector, as well as linkIDs; species are
        #aggregated by trimming links.
        res(res0|con0) = [];
        con(res0|con0) = [];
        linkTypes(res0|con0) = [];

        #My worry here is that a concomittant (or parsitic) link gets
        #aggregated with a non-concomittant (or parsitic) link.  That
        #manifests by deleting the particular superfluous link with the
        #link ID that I (may) want to keep.  Could go through the same
        #rigamarole as I do with the species names for the links, but that
        #would require another loop in here and honestly... I am too Lazy
        #to figure that out.
        #Adding concomittant links could cause more aggregation
        #(species that prey on parasites and hosts of parasites could look
        #like species that prey on hosts of parasites and concomittantly
        #prey on their parasites themselves.  But is that an important
        #distinction?)
        linkIDs(res0|con0) = [];

        #Delete the slots corresponding to superflous species.
        para(I) = [];

        #Renewing IDs - this is where the number of species gets trimmed
        oldIds = deleteVector(deleteVector>0);
        newIds = 1:length(oldIds);

        #The replacement vector works by placing the new ID at the old ID's
        #position.
        replacementVector = deleteVector;
        replacementVector(oldIds) = newIds;

        res = replacementVector(res);
        con = replacementVector(con);
        linkPara = (para(con))&(~para(res));

        #We also need to update the species names and body sizes.
        #This is not as straightforward as above, since all equivalent
        #species need to be combined in some way.

        #Combining the names
        #Averaging the body sizes.  Weight is the total biomass of each
        #species.
        I0 = I;
        J0 = J;
        for kk = I'
            #Equivalence class of the first species in array I:
            equivClasskk = [kk;J(I==kk)];


            #now, average the bodySizes.
            equivBodySizes = bodySizes(equivClasskk);
            okEquivBodySizes = isfinite(equivBodySizes)&(equivBodySizes>0);
            geo_means = geomean(equivBodySizes(okEquivBodySizes),2);
            if isempty(geo_means)
                geo_means = nan;
            end
            bodySizes(equivClasskk(end)) = geo_means;
            biomasses(equivClasskk(end)) = nansum(biomasses(equivClasskk));

            bodySizes(equivClasskk(1:end-1)) = 0;
            biomasses(equivClasskk(1:end-1)) = 0;
            for ll = equivClasskk(1:end-1)'
                #Need the if statement in case equivClasskk is empty.
                if numel(ll)>0
                    #Saving names
                    newNamesCell{ii}(equivClasskk(end)) = ...
                        strcat(newNamesCell{ii}(equivClasskk(end)),...
                        ';',newNamesCell{ii}(ll));
                    newClassCell{ii}(equivClasskk(end)) = ...
                        strcat(newClassCell{ii}(equivClasskk(end)),...
                        ';',newClassCell{ii}(ll));
                    J(I==ll) = [];
                    I(I==ll) = [];
                end
            end


        end

        bodySizes(I0) = [];
        biomasses(I0) = [];
        newNamesCell{ii}(I0) = [];
        newClassCell{ii}(I0) = [];
    end
    
    totSpecies(ii) = length(newClassCell{ii});
    MxTypes = sparse(res,con,linkTypes,totSpecies(ii),totSpecies(ii));
    [res,con,linkTypes] = find(MxTypes);
    linkPara = para(con)&(~para(res));
    linkParaCell{ii} = linkPara;
    linkListCell{ii} = [res con];
    
    totPara(ii) = sum(para);
    connectances(ii) = length(res)/totSpecies(ii)^2;
    Ls(ii) = length(res);

    properties = zeros(length(bodySizes),12);
    properties(:,1:11) = calculateLocalProperties(res,con);
    properties(:,12) = bodySizes;
    properties(:,13) = biomasses;
    propertiesCell{ii,1} = properties;
    propertiesCell{ii,2} = para;

    binParaCell{ii} = para;
    binEndoCell{ii} = zeros(totSpecies(ii),1);
    binEctoVertCell{ii} = zeros(totSpecies(ii),1);
    binEctoInvertCell{ii} = zeros(totSpecies(ii),1);
    binAutoCell{ii} = zeros(totSpecies(ii),1);


    #Ectotherm Vertebrates: one of these classes
    binAgna = ~cellfun('isempty',strfind(newClassCell{ii},'Agnatha'));
    binChon = ~cellfun('isempty',strfind(newClassCell{ii},'Chondrichthyes'));
    binOsth = ~cellfun('isempty',strfind(newClassCell{ii},'Ostheichthyes'));
    binOste = ~cellfun('isempty',strfind(newClassCell{ii},'Osteichthyes'));
    binAmph = ~cellfun('isempty',strfind(newClassCell{ii},'Amphibia'));
    binRept = ~cellfun('isempty',strfind(newClassCell{ii},'Reptilia'));
    binEctoVertCell{ii} = binAgna|binChon|binOste|binAmph|binRept|binOsth;

    #Endotherms: one of these classes
    binAvesCell = ~cellfun('isempty',strfind(newClassCell{ii},'Aves'));
    binMammCell = ~cellfun('isempty',strfind(newClassCell{ii},'Mammalia'));
    binEndoCell{ii} = binAvesCell|binMammCell;
    #sum(binEndoCell{ii}&binEctoVertCell{ii})
    #Autotrophs: generality = 0
    binAutoCell{ii} = properties(:,2)==0;
    totBasal(ii) = sum(binAutoCell{ii});
    totFreeCon(ii) = totSpecies(ii) -totBasal(ii) - totPara(ii);

    #sum(binEndoCell{ii}&binAutoCell{ii})
    #sum(binEctoVertCell{ii}&binAutoCell{ii})
    #Ectotherm invertebrates: the rest.  There are tiny little critters
    #that maybe shouldn't be classified as ectotherm invertebrates (i.e.
    #single cell organisms).. there aren't very many nodes, but might they
    #be important enough to delineate?
    binEctoInvertCell{ii} = ...
        ~(binEctoVertCell{ii}|binEndoCell{ii}|binAutoCell{ii});

    speciesTypeCell{ii} = nan(totSpecies(ii),1);
    speciesTypeCell{ii}(binEctoInvertCell{ii}) = 1;
    speciesTypeCell{ii}(binEctoVertCell{ii}) = 2;
    speciesTypeCell{ii}(binEndoCell{ii}) = 3;
    speciesTypeCell{ii}(binAutoCell{ii}) = 0;

    linkTypesCell{ii} = linkTypes;
    
    if concomittantWebs
        cd('../../Data/Processed')
        fid = fopen(sprintf('%sWebCon.csv',shortWebNames{ii}),'w');
        fprintf(fid,'res,con,para,con\n');
        fprintf(fid,'%u,%u,%u,%u\n',[res,con,linkTypes]');
        fclose(fid);
        cd(codeDir)

        fid = fopen(sprintf('%sPropertiesCon.csv',shortWebNames{ii}),'w');
        #                1         2    3         4             5           6               7           8               9           10      11
        #varNames = {'clustCoef','gen','vul','meanVulPrey','meanImpPrey','meanGenPred','meanImpPred','minSPToBasal','numConnBasal','SWTL','inLoop'};
        fprintf(fid,'\n');
        fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[properties para]');
        fclose(fid);

    elseif topoConcomittant
        continue
    else
        cd('data/Processed')
        fid = fopen(sprintf('%sWeb.csv',shortWebNames{ii}),'w');
        fprintf(fid,'res,con,para\n');
        fprintf(fid,'%u,%u,%u\n',[res,con,linkPara]');
        fclose(fid);
        cd(codeDir)

        fid = fopen(sprintf('%sProperties.csv',shortWebNames{ii}),'w');
        %                1         2    3         4             5           6               7           8               9           10      11
        %varNames = {'clustCoef','gen','vul','meanVulPrey','meanImpPrey','meanGenPred','meanImpPred','minSPToBasal','numConnBasal','SWTL','inLoop'};
        fprintf(fid,'\n');
        fprintf(fid,'%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n',[properties para]');
        fclose(fid);
    end
    
    
end

cd('data/Processed')
close all
if ectoOnly
save('webGenerationEctoInverts.mat')
else
save Generation.mat -mat7-binary
end

cd(codeDir)



