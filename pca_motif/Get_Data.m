function Get_Data(data_folder)

    Datasets = {'ALL' 'BP' 'NBN' 'NSC' 'NEUROGENESISuBPuNBN'};
	Ski = {'WT' 'KO'};
    
    for ct = Datasets
		for ski = Ski
			disp([ct{:} ' ' ski{:}])
			clear T
			
			ismara_folder = [data_folder '/' ct{:} '_' ski{:}];

			A = load([ismara_folder '/activities']);
			MeanA = load([ismara_folder '/averaged_activities']);
			DeltaMeanA = load([ismara_folder '/averaged_deltas']);

			T.All = A';
			T.Mean = MeanA';
			T.Std = DeltaMeanA;

			[~,Motifs] = textread([ismara_folder '/mat_indices'],'%u\t%s\n');
			
			samples.All = textread([ismara_folder '/sample_indices'],'%s\n');
			samples.Mean = textread([ismara_folder '/averaged_samples'],'%s\n');
			samples.Std = samples.Mean;

			save(['data/my_data_' ct{:} '_' ski{:}])
		end
    end
