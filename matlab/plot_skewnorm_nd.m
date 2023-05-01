function plot_skewnorm_nd(fig, S, us, sigmas, lambdas, ws)
%     ndims = size(S,1);
%     assert(ndims == size(us,1));
%     assert(ndims == size(sigmas,1));
%     assert(ndims == size(lambdas,1));
%     
% 	x = min(min(S)):0.01:max(max(S));
%     M = numel(w);
%     assert(M == size(us,2));
%     assert(M == size(sigmas,2));
%     assert(M == size(lambdas,2));
%     
%     for j = 1:M
%         for p = 1:ndims
%             y(p,:,j) = ws(j) * skew_norm_pdf(x, u(j), sigma(j), lambda(j));
%         end
%     end
% 	plot(ax, x, y*reshape(w, [M,1]))
end
