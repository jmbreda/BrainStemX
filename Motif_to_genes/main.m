clear all; close all; clc;

addpath('~/Matlab_Scripts/')

n_comp = 2;
lw = 2; 

% load gene expression (<log(tpm)>)
load('/scicore/home/zavolan/breda/BrainStemX/data/Matlab_Tables/Gene_SampleName_MeanLogTPM_Pseudocount_CST_gt2.mat');
Bulk = MeanLogBulk;
clear MeanLogBulk;
%load('/scicore/home/zavolan/breda/BrainStemX/data/Matlab_Tables/Gene_SampleName_StdDevLogTPM_Pseudocount_CST_gt2.mat');

% Keep only WT
idx = find(~cellfun(@isempty, regexp(Bulk.Properties.VariableNames,'KO') ));
Bulk(:,idx) = [];

% Get Gene name
tmp = cellfun(@(g)strsplit(g,'|'),Bulk.GeneID,'UniformOutput',0);
tmp = [tmp{:}];
Bulk.GeneID = tmp(2:2:end)';

% Accumulate Rows with same gene name
[GeneUniq,~,ind_uniq] = unique(Bulk.GeneID);
[xx, yy] = ndgrid(ind_uniq,1:size(Bulk,2)-1);
M = Bulk{:,2:end};
M = accumarray([xx(:) yy(:)],M(:));
var_name = Bulk.Properties.VariableNames(2:end);
Bulk = [table(GeneUniq) array2table(M)];
Bulk.Properties.VariableNames = ['GeneID' var_name];


celltype = {'NSC','BP','NBN','NEUROGENESISuBPuNBN'};
for ct = celltype
    disp(ct{:})
    disp(n_comp)


    % Keep only celltype
    switch ct{:}
    case 'NEUROGENESISuBPuNBN'
        idx = [4:8 12:24];
    otherwise
        idx = find(~cellfun(@isempty, regexp(Bulk.Properties.VariableNames,ct{:}) ));
    end
    E = Bulk{:,idx};
    E_genes = Bulk.GeneID;

    % load motif activity
    ismara = ['/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna/' ct{:} '_WT'];
    A = load([ismara '/averaged_activities']);
    Adj = load([ismara '/Adjacency_motif_gene']);
    [~,motifs] = textread([ismara '/mat_indices'],'%u\t%s\n');
    genes = textread([ismara '/gene_indices'],'%s\n');

    [coeff_e,score_e,latent_e] = pca(E');
    f_v_g = coeff_e(:,1:n_comp).^2*latent_e(1:n_comp)/sum(latent_e);
    [f_v_g,idx] = sort(f_v_g,'descend');
    E_genes = E_genes(idx);

    [coeff_a,score_a,latent_a] = pca(A);
    f_v_m = coeff_a(:,1:n_comp).^2*latent_a(1:n_comp)/sum(latent_a);
    [f_v_m,idx] = sort(f_v_m,'descend');
    motifs = motifs(idx);
    Adj = Adj(idx,:);
    A = A(:,idx);

    figure
    plot(linspace(0,1,length(f_v_g)),cumsum(f_v_g),'Linewidth',lw)
    hold on 
    plot(linspace(0,1,length(f_v_m)),cumsum(f_v_m),'Linewidth',lw)
    
    gene_mot_score = sum(exp(Adj).*repmat(f_v_m,1,size(Adj,2)))'; 
    [gene_mot_score,idx] = sort(gene_mot_score,'descend');
    genes = genes(idx);

    X = 10:10:15000;
    inter_g = zeros(size(X));
    i = 0;
    %inter_m = zeros(size(X));
    for k=X
        if ~mod(k,1e3)
            disp(k)
        end
        i=i+1;
        inter_g(i) = length(intersect(genes(1:k),E_genes(1:k)))/k;
    end

    plot(X/max(X),inter_g,'Linewidth',lw)
    xlabel('genes/motifs frac. sorted by variance')
    dim = [8 6];
    set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
    print(gcf,['Fig/Gene_overlap_' num2str(n_comp) '_' strrep(ct{:},'|','_')],'-dpdf');
end


figHandle = figure;
p1 = plot([1:10], [1:10]);
hold on;
p2 = plot([1:10], [1:10]);
p3 = plot([1:10], [1:10]);

legHandle = legend('fraction of total variance in gene','fraction of total variance in motif',...
'overlap fraction genes and motifs');

saveLegendToImage(figHandle, legHandle, 'Fig/gene_overlap_legend', 'pdf');

