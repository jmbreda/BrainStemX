
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

% Compute pca on WT only :
idx_WT = find(~cellfun( @isempty,regexp(samples.Mean,'WT') ));
[coeff,score_mean,latent,tsquared] = pca(T.Mean(:,idx_WT)');
if inv_pc1
	coeff(:,1) = -1*coeff(:,1);
	score_mean(:,1) = -1*score_mean(:,1);
end
if inv_pc2
    coeff(:,2) = -1*coeff(:,2);
    score_mean(:,2) = -1*score_mean(:,2);
end

% Project motif activity of replicates on the replicate average pc :
score = [T.All - repmat(mean(T.All,2),1,length(samples.All))]' * coeff(:,1:2);

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

%my_axis(4) = axis(1.1*my_axis)
xlabel(['pc_1 : ' num2str(latent(1)/sum(latent),2)]);
ylabel(['pc_2 : ' num2str(latent(2)/sum(latent),2)]);
%set(gca,'Xtick',[],'Ytick',[])
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0  0 dim],'PaperSize',[dim]);
print(gcf,['Fig/motif_pca_' celltype],'-dpdf');
%print(gcf,['Fig/motif_pca_' celltype '_maria'],'-dpdf');

