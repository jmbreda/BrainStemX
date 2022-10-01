
if strcmp(celltype,'ALL')
    load('data/my_data_NSC.mat')
    NSC = T.Mean;
    [coeff,~,~,~] = pca(NSC');
    t_nsc = coeff(:,1:2);
end

load(['data/my_data_' celltype '.mat'])

my_color = [];
timepoint = {};
for ii = fieldnames(T)'
    for i = 1:length(samples.(ii{:}))
        tmp = strsplit(samples.(ii{:}){i},'_');
        my_color.(ii{:})(i,:) = C( find(strcmp(tmp{1},Timepoint)) ,:);
        timepoint{end+1} = tmp{1};
    end
end
timepoint= unique(timepoint);
Nt = length(timepoint);

if strcmp(celltype,'ALL')
    T_Mean_centered = T.Mean - repmat(mean(T.Mean,2),1,length(samples.Mean));
    for d = 1:size(T_Mean_centered,2);
        % substract the projection of M on t
        for i = 1:size(t_nsc,2)
            T_Mean_centered(:,d) = T_Mean_centered(:,d) - (T_Mean_centered(:,d)'*t_nsc(:,i))*t_nsc(:,i);
        end
    end
    T.Mean = T_Mean_centered;
end

[coeff,score_mean,latent,tsquared] = pca(T.Mean');
if inv_pc1
	coeff(:,1) = -1*coeff(:,1);
	score_mean(:,1) = -1*score_mean(:,1);
end
if inv_pc2
    coeff(:,2) = -1*coeff(:,2);
    score_mean(:,2) = -1*score_mean(:,2);
end
[~,ind_proj] = sort(coeff(:,1:2).^2 * latent(1:2),'descend');
out_motifs = {'Hdx','Nfatc1','Irf2_Irf1_Irf8_Irf9_Irf7','Sp100'};
[~,idx_out,~] = intersect(Motifs,out_motifs,'stable');
ind_proj = setdiff(ind_proj,idx_out,'stable');

top_coeff = coeff(ind_proj(1:n_best),1:2);
tree = linkage(top_coeff,'ward');
my_coeff_ind = cluster(tree,'maxclust',n_clust);
for ii = 1:n_clust
	my_coeff(ii,:) = mean(top_coeff(my_coeff_ind==ii,:),1);
end

% Project motif activity of replicates on the replicate average pc :
score = [T.All - repmat(mean(T.All,2),1,length(samples.All))]' * coeff(:,1:2);



figure('visible','off')
hold on;
for k = unique(my_coeff_ind)'
    plot([0 l_proj*my_coeff(k,1)],[0 l_proj*my_coeff(k,2)],'-','color',my_grey)
    h = text(l_proj*my_coeff(k,1),l_proj*my_coeff(k,2),Motifs(ind_proj(my_coeff_ind==k)),'Fontsize',default_fs,'interpreter','none');
    if my_coeff(k,1) < 0
        set(h,'HorizontalAlignment','right');
    end
    tmp = Motifs{ind_proj(my_coeff_ind==k)};
    if strcmp(tmp,'Pou2f2_Pou3f1') && strcmp(celltype,'NEUROGENESISuBPuNBN')
        set(h,'VerticalAlignment','top');
    end
    if strcmp(tmp,'E2f1') && strcmp(celltype,'NEUROGENESISuBPuNBN')
        set(h,'VerticalAlignment','bottom');
    end
    if strcmp(tmp,'E2f2_E2f5') && strcmp(celltype,'NEUROGENESISuBPuNBN')
        set(h,'VerticalAlignment','top');
    end

    if strcmp(tmp,'Sox2') && strcmp(celltype,'NSC')
        set(h,'VerticalAlignment','top');
    end
    if strcmp(tmp,'Stat2') && strcmp(celltype,'NSC')
        set(h,'VerticalAlignment','bottom');
    end

	if strcmp(celltype,'NSC_WT')
		switch tmp
		case 'Stat2'
			set(h,'VerticalAlignment','bottom');
		end
	end

	if strcmp(celltype,'NSC_KO')
		switch tmp
		case 'Tbx4'
			set(h,'VerticalAlignment','bottom');
		end
	end

	if strcmp(celltype,'BP_KO')
		switch tmp
		case 'Foxd1'
			set(h,'VerticalAlignment','bottom','HorizontalAlignment','center');
		end
	end

	if strcmp(celltype,'NBN_KO')
		switch tmp
		case 'Spib'
			set(h,'VerticalAlignment','bottom');
		end
	end
end
axis(my_axis)
xlabel(['pc_1 : ' num2str(latent(1)/sum(latent),2)]);
ylabel(['pc_2 : ' num2str(latent(2)/sum(latent),2)]);
set(gca,'Xtick',[],'Ytick',[])
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0  0 dim],'PaperSize',[dim]);
print(gcf,['Fig/motif_projection_' celltype],'-dpdf');
%print(gcf,['Fig/motif_projection_' celltype '_maria'],'-dpdf');





figure('visible','off')
hold on
my_axis = [0 0 0 0];
for ct = 1:length(CellTypes)
    idx = find(~cellfun(@isempty,regexp(samples.All,CellTypes{ct})));
    if ~isempty(idx)
        my_marker = Markers{ct};
        scatter(score(idx,1),score(idx,2),30,my_color.All(idx,:),my_marker,'filled')

		my_axis(1) = min(my_axis(1),min(score(idx,1)));
        my_axis(2) = max(my_axis(2),max(score(idx,1)));
        my_axis(3) = min(my_axis(3),min(score(idx,2)));
        my_axis(4) = max(my_axis(4),max(score(idx,2)));
    end
end
if strcmp(celltype,'NSC')
    my_xlim = get(gca,'xlim');
    my_ylim = get(gca,'xlim');

    idx.E = 1:2;
    idx.N = 3:7;
    idx.G = 8:10;
    
    for p = {'E','N','G'}
        Score.(p{:}) = mean(score_mean(idx.(p{:}),1:2),1);
        Score.(p{:}) = Score.(p{:})/norm(Score.(p{:}));
        phi.(p{:}) = atan2( Score.(p{:})(1), Score.(p{:})(2) );
    end
    plot(0,0,'.','MarkerSize',1)
end
axis(1.1*my_axis)
xlabel(['pc_1 : ' num2str(latent(1)/sum(latent),2)]);
ylabel(['pc_2 : ' num2str(latent(2)/sum(latent),2)]);
set(gca,'Xtick',[],'Ytick',[])
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0  0 dim],'PaperSize',[dim]);
print(gcf,['Fig/motif_pca_' celltype],'-dpdf');
%print(gcf,['Fig/motif_pca_' celltype '_maria'],'-dpdf');

fid = fopen(['output/motifs_score_1_2_' celltype '.txt'],'w');
for i = 1:length(samples.All)
    fprintf(fid,'%s\t%f\t%f\n',samples.All{i},score(i,1),score(i,2));
end
fclose(fid);

var_pc_1_2 = coeff(:,1:2).^2 * (latent(1:2)/sum(latent));
fid = fopen(['output/motifs_var_' celltype '.txt'],'w');
for i = ind_proj'
    fprintf(fid,'%s\t%f\n',Motifs{i},var_pc_1_2(i));
end
fclose(fid);

