function plot_skewnorm(x, alpha, beta, w)
	x = [min(x):0.01:max(x)];
	y = gampdf(x, alpha, 1/beta);
	plot(x, y*w)
end
