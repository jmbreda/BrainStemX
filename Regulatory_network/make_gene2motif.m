clear all;

ismara_folder = '/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna/ALL_ALL';
[~,motifs] = textread([ismara_folder '/mat_indices'],'%u\t%s\n');

% Get TF genes
k = 1;
for m = motifs'
    for g = strsplit(m{:},'_')
        my_genes{k,1} = g{:};
        my_genes_motif{k,1} = m{:};
        k = k+1;
    end
end
gene2motif = containers.Map(my_genes,my_genes_motif);
save('gene2motif.mat','gene2motif');
