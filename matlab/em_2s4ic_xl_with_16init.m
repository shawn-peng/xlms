
loaddata_xl;

mat2 = mat(:,mat(2,:)~=0);
best_ll = -inf;

s1 = mat2(1,:);
s1 = sort(s1, 'descend');
q1 = quantile(s1, 0.01);
s1 = s1(:, (s1 > q1));

info_file = sprintf('%s%s.json', '../results/info/', dataset_name);
info = load_json(info_file)

% d = sprintf('../figures/%s', dataset_name);
% if ~exist(d)
%     mkdir(d);
% else
%     delete([d, '/*']);
% end

for j = 0:16
    sl1 = get_sign(j, 1);
    sl2 = get_sign(j, 2);
    sl3 = get_sign(j, 3);
    sl4 = get_sign(j, 4);
    if sl1 < 0
        continue
    end
    try
        [params, ll, ll1] = EM2_4ic_xl(mat2,sl1,sl2,sl3,sl4,'tolerance',1e-4);
    catch ME
        continue
    end
    if ll > best_ll
        best_ll = ll;
        best_params = params;
        best_sls = [sl1, sl2, sl3, sl4];
    end
    
    close();
end


ws = best_params{1};
theta = best_params{2};
[u_c, sigma_c, lambda_c] = unpack_skntheta(theta(1));
[u_ic, sigma_ic, lambda_ic] = unpack_skntheta(theta(2));
[u_i1, sigma_i1, lambda_i1] = unpack_skntheta(theta(3));
[u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta(4));

sl1 = best_sls(1);
sl2 = best_sls(2);
sl3 = best_sls(3);
sl4 = best_sls(4);


function s = get_sign(j, d)
if j == 0
    s = 0;
elseif bitand(j-1, bitshift(1, d-1))
    s = -1;
else
    s = 1;
end
end

