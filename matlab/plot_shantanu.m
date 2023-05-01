
list_species = {
%     'HeLa01ng'
    'HeLa1ng'
%     'HeLa10ng'
%     'HeLa50ng'
    'HeLa100ng'
};

list_methods = {
    'SNMax1'
};


n_sp = size(list_species, 1);

n_m = size(list_methods, 1);


for sp_i = 1:n_sp
    species = list_species{sp_i}
    for me_j = 1:n_m
        method = list_methods{me_j};
        plot_dist_obj
    end
end
        