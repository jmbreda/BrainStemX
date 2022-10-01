clear all; close all; clc;

ismara_folder = '/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna/ALL_WT';

E = load([ismara_folder '/prom_expression']);
N = load([ismara_folder '/sitecount_mat']);
A = load([ismara_folder '/activities']);
[~,Prom] = textread([ismara_folder '/prom_indices'],'%u\t%s\n');

my_genes = textread('my_genes.txt','%s');
my_prom = {};
gene2prom = '/scicore/home/nimwegen/GROUP/software/mara_data/mm10_f5/clusters_150_gene_associations';
for g = my_genes'
    [stat,out] = system(['grep ' g{:} ' ' gene2prom ' | cut -f6']);
    my_prom{end+1} = out(1:end-1);
end

[~,idx_prom,~] = intersect(Prom,my_prom,'stable');
E = E(idx_prom,:);
N = N(idx_prom,:);

[N_prom,N_samp] = size(E);
[~,N_mot] = size(N);

%Normalizing
N = N - repmat(mean(N,1),N_prom,1);
A = A - repmat(mean(A,1),N_samp,1);
E = E - repmat(mean(E,1),N_prom,1) - repmat(mean(E,2),1,N_samp) + repmat(mean(mean(E)),N_prom,N_samp);

E_fit = N*A';

C = zeros(N_prom,1);
for i = 1:N_prom
    C(i) = corr(E(i,:)',E_fit(i,:)');
end


figure;
plot(E',E_fit','.')
