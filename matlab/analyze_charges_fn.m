function analyze_charges_fn(species)
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