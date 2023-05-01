
list_species = {
%     'HeLa01ng'
%     'HeLa1ng'
%     'HeLa10ng'
    'HeLa50ng'
%     'HeLa100ng'
};
% list_species = {
%     'c_elegans'
%     'drosophila'
%     'e_coli'
%     'human'
%     'mouse'
% }
results_folder = 'test_search/est_results/';
% results_folder = 'test_search/est_results_nist/';

list_methods = {
    'SNMax1'
};


n_sp = size(list_species, 1);

n_m = size(list_methods, 1);

results_arr = {};

for sp_i = 1:n_sp
    species = list_species{sp_i}
    results_arr = {}
    for me_j = 1:n_m
%         run parse_psm.m
        method = list_methods{me_j};
        
        load(['test_search/shantanu/',species,method,'_B.mat']);

        obj = {};
        obj.ds = S.ds;
        obj.algo = (S.algo);
        obj.s1 = cell2mat(S.s1);
        obj.pdfM = cell2mat(S.pdfM);
        obj.pdfI1 = cell2mat(S.pdfI1);
        obj.pdfC = cell2mat(S.pdfC);
        obj.cdfTrue = cell2mat(S.cdfTrue);
        obj.cdf = cell2mat(S.cdf);
        obj.fdr = cell2mat(S.fdr);
        obj.t1p = (S.t1p);
        % obj.t01p = S.;
        % obj.t10p = thres10;
        obj.ll = (S.ll1);
        obj.deltaCdf = (S.deltaCdf);
        obj.pi_C = (S.pi_C);

        obj.t1p

        results_arr{end+1} = obj;
        
        plot_dist_obj
    end
    json = jsonencode(results_arr);
    
    fid = fopen(['test_search/shantanu/json/', species, '.json'],'wt');
    fprintf(fid, json);
    fclose(fid);
end

