load('-v7','data/Processed/Generation.mat', 'OGSpeciesLinks')

# ORIGINAL SPECIES WEBS ARE SAVED IN OGSpeciesLinks:
#  OGSpeciesLinks{:,1} are the link lists: [res, con]
#  OG SpeciesLinks{:,2} are the parasite species

# Need to do a little bit of preprocessing to get everything as we expect it:
n_webs = numel(OGSpeciesLinks);
n_cells = cell(1,n_webs);
local_props = n_cells;
global_props = n_cells;
local_mean_props = n_cells;
count = 0;

Ss = zeros(n_webs, 1);
Cs = zeros(n_webs, 1);
Sfs = zeros(n_webs, 1);
Sps = zeros(n_webs, 1);
Scs = zeros(n_webs, 1);
for web_data = OGSpeciesLinks'
  count += 1
  # First, convert the string numbers in the link list to number numbers.
  web_data = web_data{1};
  LL = web_data{1};
  paras = web_data{2};
  LL = cellfun(@(x) str2num(x), LL);
  paras = cellfun(@(x) str2num(x), paras);
  sp_ids = unique(LL);
  para_ids = unique(paras);
  reindex = containers.Map(sp_ids,1:numel(sp_ids));
  res = arrayfun(@(x) reindex(x), LL(:,1));
  con = arrayfun(@(x) reindex(x), LL(:,2));
  para_ids = arrayfun(@(x) reindex(x), para_ids);
  S = numel(sp_ids);
  para = false(S,1);
  para(para_ids) = true;
  
  min_cluster_distance = 0;
  node_sizes = ones(S,1);
  counter = 1;
  A = full(sparse(res,con,1,S,S));
  [local_props_, global_props_, local_mean_props_, res_, con_] = agglomProps(A,para,node_sizes,0,1);
  local_props{count} = local_props_';
  global_props{count} = global_props_';
  local_mean_props{count} = local_mean_props_';
  Ss(count) = S;
  Cs(count) = numel(res)/S^2;
  Sfs(count) = S - sum(para);
  Sps(count) = sum(para);
  Scs(count) = sum(local_props_(:,11));
endfor
OG_global_props = [global_props{:}]';
OG_local_means = [local_mean_props{:}]';
OG_local_props = [local_props{:}]';
save -v7 og_local_means.mat OG_local_means OG_global_props OG_local_props

