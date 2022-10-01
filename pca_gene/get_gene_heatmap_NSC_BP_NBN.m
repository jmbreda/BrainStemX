

my_color = [];
timepoint = {};
for ii = fieldnames(Tab)'
    T.(ii{:}) = Tab.(ii{:})(:,celltype_ind.(ii{:}).(celltype));
    samples.(ii{:}) = Samples.(ii{:})(celltype_ind.(ii{:}).(celltype));

    for i = 1:length(samples.(ii{:}))
        tmp = strsplit(samples.(ii{:}){i},'_');
        my_color.(ii{:})(i,:) = C( find(strcmp(tmp{1},Timepoint)) ,:);
        timepoint{end+1} = tmp{1};
    end
end
timepoint= unique(timepoint);
Nt = length(timepoint);

[coeff,score,latent,tsquared] = pca(T.All');
[~,ind_proj] = sort(coeff(:,1:2).^2 * latent(1:2),'descend');
%[my_coeff_ind,my_coeff] = kmeans(coeff(ind_proj(1:n_best),1:2),n_clust,'Replicate',10);

for ii = 1:length(samples.All)
    tmp = strsplit(samples.All{ii},'_');
    sample_name{ii} = [tmp{1} ' ' tmp{3}];
end

figure
cg = clustergram(T.All(ind_proj(1:n_genes),:),'RowLabels',Genes(ind_proj(1:n_genes)),'ColumnLabels',sample_name)
cg.Colormap = redbluecmap;
