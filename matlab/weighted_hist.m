function weighted_hist(x, weight, nbins)
	[histw, intervals] = histwc(x, weight, nbins);
	w = intervals(2) - intervals(1);
	n = prod(size(x));
	bar(intervals, histw/n/w)
end

