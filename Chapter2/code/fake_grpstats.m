function varargout = fake_grpstats(X,grps,fcns)
  
  grp_set = unique(grps);
  n_grps = numel(grp_set);
  n_stats = numel(fcns);
  varargout = cell(n_stats + 1,1);
  out_array = zeros(n_grps,size(X,2));
  varargout(:) = deal(out_array);
  fn_count = 0;
  for fn = fcns
    fn_count += 1;
    fn = fn{1};
    count = 0;
    for grp = grp_set'
      count += 1;
      this_grp = grps == grp;
      varargout{fn_count}(count,:) = fn(X(this_grp,:));
    endfor
  endfor
  varargout{end} = grp_set;
endfunction
