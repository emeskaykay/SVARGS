k = 100;
nsamples = 16;
run_parallel = true;

tmin = 150;
tmax = 450;
tstep = 50;
T = (tmin : tstep : tmax);

mincoeff = 0.2;
maxcoeff = Inf;

sample_first_index = 1;
sample_last_index = 3;

data = cell(1, nsamples);
models(nsamples) = struct('opsparse', [], 'sigma', []);
pfit = zeros(1,nsamples);

if run_parallel
	poolobj = gcp('nocreate');
	if isempty(poolobj)
		parpool; % create parallel pool
	end
end

found = false(1,nsamples);
parfor j = 1 : nsamples
	found(j) = false;
	p = 5 + round(3*rand);
	
	pfactorords = [];
	switch (p)
		case 5, pfactorords = [2,3];
		case 6, pfactorords = [1,2,3];
		case 7, pfactorords = [2,5];
		case 8, pfactorords = [3,5];
		otherwise
			assert(false);
	end
	
	ccount_1 = k*(p/2 + 2*rand);
	ccount_2 = k*(0.5 + 0.5*rand)*tmin/25;
	ccount_3 = (k^2)*p/5;
	ccount = min([ccount_1, ccount_2, ccount_3]);
	ccount = round(ccount);
	
	[models(j), found(j)] = generate_sparse_model(k, p, ccount, [mincoeff, maxcoeff], pfactorords, ...
		'runparallel', false, 'maxruntime', 10, 'nmaxtries', 20);
	
	data{j} = generate_data(models(j), tmax);
	pfit(j) = p + 2;
end

numfound = nnz(found);
disp(['Found ', num2str(numfound) ' out of ', num2str(nsamples)]);

parfor j = 1 : nsamples
	dispmasker = Display_masker();
	dispmasker.Flag('char', true);
	dispmasker.Flag('struct', true);
	for n = 1 : numel(T)
		t = T(n);
		yy = data{j}(:, 1:t);
		
		disp(['Data length: ', num2str(t)]);
		disp(['Sample no.: ', num2str(j)]);
		
		lags = 1 : pfit(j);
		parameters = Svargs_params(k, lags);
		fitted_models(j) = svargs_run(yy, parameters, 1); %#ok<SAGROW>
	end
end

summary = AssessModels(models, fitted_models);
summary
