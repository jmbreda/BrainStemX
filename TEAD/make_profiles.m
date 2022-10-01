clear all; close all; clc;

A = load('../../data/Ismara/Pseudocount_CST_gt2_no_mirna/NSC_WT/averaged_activities');
dA = load('../../data/Ismara/Pseudocount_CST_gt2_no_mirna/NSC_WT/averaged_deltas');
samp = textread('../../data/Ismara/Pseudocount_CST_gt2_no_mirna/NSC_WT/averaged_samples','%s\n');
[~,motifs] = textread('../../data/Ismara/Pseudocount_CST_gt2_no_mirna/NSC_WT/mat_indices','%u\t%s\n');

my_motifs = {'Creb3l2','Crem_Jdp2','Nfkb2','Nr2f1_Nr4a1','Pitx2_Otx2','Rest','Sox6_Sox9','Tead1','Tead3_Tead4'};
for m = my_motifs
    idx_m = find(strcmp(motifs,m))

    figure('visible','off')
    errorbar(1:10,A(:,idx_m),dA(:,idx_m),'k.-')
    xlim([.5 10.5])
    set(gca,'Xtick',[],'Ytick',[]);
    dim = [4 3];
    title(m{:},'interpreter','none');
    set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
    'PaperPosition',[0 0 dim],'PaperSize',[dim])
    print(gcf,['Fig/' m{:}],'-dsvg');
end

load('../../data/Matlab_Tables/Gene_SampleName_MeanLogTPM_Pseudocount_CST_gt0.mat');
load('../../data/Matlab_Tables/Gene_SampleName_StdDevLogTPM_Pseudocount_CST_gt0.mat');

[~,idx_out] = setdiff(MeanLogBulk.Properties.VariableNames,samp);
idx_out(find(idx_out==1))= [];
MeanLogBulk(:,idx_out) = [];
StdLogBulk(:,idx_out) = [];

idx_g = find(~cellfun(@isempty,regexp(MeanLogBulk.GeneID,'Tead')));

for g = idx_g'
    disp(g)
    tmp = strsplit(MeanLogBulk.GeneID{g},'|');

    figure('visible','off')
    errorbar(1:10,MeanLogBulk{g,2:end},StdLogBulk{g,2:end},'k.-')
    xlim([.5 10.5])
    %axis off
    %set(gca,'Xtick',[],'Ytick',[]);
    title(tmp{2},'interpreter','none');
    dim = [4 3];
    set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto',...
    'PaperPosition',[0 0 dim],'PaperSize',[dim])
    %print(gcf,['Fig/' tmp{2}],'-dsvg');
    print(gcf,['Fig/' tmp{2}],'-dpdf');
end

