function analyze_results_fn(species)
	species

	% data_dir = 'test_search/matdata_nist/';
	% results_folder = 'test_search/est_results_nist_allinits/';

	% data_dir = 'test_search/matdata_hela/';
	% results_folder = 'test_search/est_results_hela_allinits/';

	data_dir = 'test_search/matdata_pride/';
	results_folder = 'test_search/est_results_pride_allinits/';

	xlim_high = 65;
	plotcdf = true;
	plotthres = false;

	figw = 240;
	figh = 120;

	c_color = [0.8500 0.3250 0.0980];
	i1_color = [0.9290 0.6940 0.1250];
	i2_color = [0 0.4470 0.7410];
	mix_color = [0.4940 0.1840 0.5560];

	cdfcolor = [0.6350 0.0780 0.1840];
	legending = false;
	legending2 = false;

	linewidth = 1;

    load([data_dir,species,'_data.mat'])
    run_all
end
