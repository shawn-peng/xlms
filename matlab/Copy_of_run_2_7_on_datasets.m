
datasets = {
    'ecoli_xl';
    'MS2000225';
    'alban';
    'RPA';
    'CPSF';
    'Alinden';
    'ALott';
    'KKT4';
    'QE';
    'peplib';
};

for i = 1:size(datasets)
    dataset_name = datasets{i};
    load(['../results/matdata/scoremats/', dataset_name, '.mat']);
    EM2_7j_xl(mat,1,-1,-1,-1)
    title(dataset_name);
end