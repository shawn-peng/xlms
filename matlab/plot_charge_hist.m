
species_o = species;

species_folder = [results_dir, species, '/'];
datadir = 'test_search/matdata_hela/';

figure
hold on;

legends = {};
for c = 2:3
    species = species_o;
    datafile = [datadir, species, '.c', num2str(c), '_data.mat'];
    load(datafile);
    s1 = omat(:,1);
    s1 = s1(s1~=0);

    histogram(s1, 'Normalization', 'pdf');
    xlabel('-log(EValue)')
    ylabel('PDF')
    
    legends{end+1} = ['Charge ', num2str(c)];

end
legend(legends);
species = species_o;

saveas(gcf, [species_folder, 'charges_hist.png']);

