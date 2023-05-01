function score_thres = fdr_thres(s, fdrs, thres)

assert(all(size(s) == size(fdrs)));

[s, ind] = sort(s,'descend');
fdrs = fdrs(ind);
thres_ind = find(fdrs >= thres, 1);
score_thres = s(thres_ind);

end