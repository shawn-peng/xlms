
load([data_dir,species,'_data.mat'])
species_folder = [results_folder,species];

for me_j = 1:n_m
    method = list_methods{me_j};

    load([species_folder, '/params/', method, '.mat'])

    if strcmp(method, '_1s2c')
        algo = '1S2D skew normal';
        alpha = theta.alpha;
        u_c = theta.theta_c.u;
        sigma_c = theta.theta_c.sigma;
        lambda_c = theta.theta_c.lambda;
        u_i = theta.theta_i.u;
        sigma_i = theta.theta_i.sigma;
        lambda_i = theta.theta_i.lambda;
    elseif strcmp(method, '_1s2ca')
        algo = '1S2D gamma & gaussian';
        alpha = theta.alpha;
        u_c = theta.theta_c.u;
        sigma_c = theta.theta_c.sigma;
        a_i = theta.theta_i.a;
        b_i = theta.theta_i.b;
        gamma_i = theta.theta_i.gamma;
    elseif strcmp(method, '_2s3ci') || strcmp(method, '_2s3ct')
        algo = '2S3D skew normal';
        if strcmp(method, '_2s3ct')
            algo = [algo, ' truncated'];
        end
        alpha = theta.alpha;
        beta = theta.beta;
        [u_c, sigma_c, lambda_c] = unpack_skntheta(theta.theta_c);
        [u_i, sigma_i, lambda_i] = unpack_skntheta(theta.theta_i);
        [u_i2, sigma_i2, lambda_i2] = unpack_skntheta(theta.theta_i2);
    end

    run(['result_curv', method])

    obj = {};
    obj.ds = species;
    obj.algo = algo;
    obj.s1 = s1;
    obj.pdfM = y1;
    obj.pdfI1 = yi1;
    obj.pdfC = yc;
    obj.cdfTrue = h1emp;
    obj.cdf = h1;
    size(obj.cdf)
    obj.fdr = fdr_curv;
    obj.t1p = thres;
    obj.t01p = thres01;
    obj.t10p = thres10;
    obj.ll = ll_1;
    obj.deltaCdf = sdcdf;
    obj.pi_C = alpha;

    obj.t1p
    figure
    hold on;
    histogram(s1, 'Normalization', 'pdf')
    plot(s1, y1)
    plot(s1, yi1)
    plot(s1, yc)
    plot(s1, h1emp)
    plot(s1, h1)
    plot(s1, fdr_curv)
    results_arr{end+1} = obj;
end
json = jsonencode(results_arr);

fid = fopen([json_dir, species, '.json'],'wt');
fprintf(fid, json);
fclose(fid);
results_arr = {};
