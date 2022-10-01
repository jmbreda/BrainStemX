
NSC = Tab.All(:,celltype_ind.All.NSC);
[coeff,~,~,~] = pca(NSC');
t_nsc = coeff(:,1:2);    

my_color = [];
timepoint = {};
for ii = fieldnames(Tab)'
    ind = find( ~cellfun(@isempty,regexp(Samples.(ii{:}),'WT')) );
    T.(ii{:}) = Tab.(ii{:})(:,ind);
    samples.(ii{:}) = Samples.(ii{:})(ind);

    for i = 1:length(samples.(ii{:}))
        tmp = strsplit(samples.(ii{:}){i},'_');
        my_color.(ii{:})(i,:) = C( find(strcmp(tmp{1},Timepoint)) ,:);
        timepoint{end+1} = tmp{1};
    end
end 
timepoint = unique(timepoint);
Nt = length(timepoint);

T_All_centered = T.All - repmat(mean(T.All,2),1,length(samples.All));
for d = 1:size(T_All_centered,2);
    % substract the projection of M on t
    for i = 1:size(t_nsc,2)
        T_All_centered(:,d) = T_All_centered(:,d) - (T_All_centered(:,d)'*t_nsc(:,i))*t_nsc(:,i);
    end
end

[coeff,score,latent,tsquared] = pca(T_All_centered');
if inv_pc1
    coeff(:,1) = -1*coeff(:,1);
    score(:,1) = -1*score(:,1);
end
[~,ind_proj] = sort(coeff(:,1:2).^2 * latent(1:2),'descend');
[my_coeff_ind,my_coeff] = kmeans(coeff(ind_proj(1:n_best),1:2),n_clust,'Replicate',10);




figure('visible','off')
hold on
for k = unique(my_coeff_ind)'
    plot([0 l_proj*my_coeff(k,1)],[0 l_proj*my_coeff(k,2)],'-','color',my_grey)
    h = text(l_proj*my_coeff(k,1),l_proj*my_coeff(k,2),Genes(ind_proj(my_coeff_ind==k)),'Fontsize',default_fs);
    if my_coeff(k,1) < 0
        set(h,'HorizontalAlignment','right');
    end
    tmp = Genes{ind_proj(my_coeff_ind==k)}
end
axis(my_axis)
xlabel(['pc_1 : ' num2str(latent(1)/sum(latent),2)]);
ylabel(['pc_2 : ' num2str(latent(2)/sum(latent),2)]);
set(gca,'Xtick',[],'Ytick',[])
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0  0 dim],        'PaperSize',[dim]);
print(gcf,['Fig/gene_projection_' celltype],'-dpdf');




figure('visible','off')
hold on
for ct = 1:length(CellTypes)
    idx = find(~cellfun(@isempty,regexp(samples.All,CellTypes{ct})));
    if ~isempty(idx)
        my_marker = Markers{ct};
        scatter(score(idx,1),score(idx,2),30,my_color.All(idx,:),my_marker,'filled')
    end
end
xlabel(['pc_1 : ' num2str(latent(1)/sum(latent),2)]);
ylabel(['pc_2 : ' num2str(latent(2)/sum(latent),2)]);
set(gca,'Xtick',[],'Ytick',[])
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0  0 dim],        'PaperSize',[dim]);
print(gcf,['Fig/gene_pca_' celltype],'-dpdf');


