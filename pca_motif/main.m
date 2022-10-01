clc; close all; clear all;

addpath('~/Matlab_Scripts/cbrewer/')

default_fs = 6;
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',default_fs)
dim = [8 6];
my_grey = [.8 .8 .8];
n_best = 20;

Timepoint = {'E105','E115','E125','E135','E145','E155','E165','E175','E185','PN'};
%C = cbrewer('div','Spectral',10); 
C = cbrewer('div','RdYlGn',10);

CellTypes = {'NSC','BP','NBN'};
Markers = {'s','o','^'};

if 0
    %dataset = '/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna_maria';
    dataset = '/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna';
    Get_Data(dataset)
end

if 0
% PCA ALL orth to NSC
celltype = 'ALL_WT';
l_proj = 3.5;
n_clust = 5;
my_axis = [-.7 .9 -.35 .75];
inv_pc1=0;
inv_pc2 = 0;
get_pca_mean_activity_mot_proj
end

if 0
celltype = 'NEUROGENESISuBPuNBN';
n_clust = 6;
l_proj = 1.2;
my_axis = [-.4 .25 -.2 .35];
inv_pc1 = 1;
inv_pc2 = 0;
get_pca_mean_activity_mot_proj
end

if 0
celltype = 'NSC_WT';
inv_pc1 = 0;
n_clust = 6;
l_proj = 3.5;
my_axis = [-.9 1 -1.2 1];
inv_pc1 = 0;
inv_pc2 = 1;
get_pca_mean_activity_mot_proj
end

if 0
celltype = 'NSC_KO';
inv_pc1 = 0;
n_clust = 7;
l_proj = 3.5;
my_axis = [-1 1 -1.2 1];
inv_pc1 = 0;
inv_pc2 = 1;
get_pca_mean_activity_mot_proj
end


if 0
celltype = 'BP_WT';
n_clust = 5;
l_proj = 3;
my_axis = [-1.1 .9 -.55 .5];
inv_pc1 = 0;
inv_pc2 = 0;
get_pca_mean_activity_mot_proj
end

if 0
celltype = 'BP_KO';
n_clust = 5;
l_proj = 3;
my_axis = [-.9 .9 -.4 .6];
inv_pc1 = 0;
inv_pc2 = 0;
get_pca_mean_activity_mot_proj
end

if 1
celltype = 'NBN_WT';
n_clust = 7;
l_proj = 3;
my_axis = [-.9 .9 -.5 .5];
inv_pc1 = 1;
inv_pc2 = 0;
figure('visible','off')
get_pca_mean_activity_mot_proj
end

if 0
celltype = 'NBN_KO';
n_clust = 3;
l_proj = 3;
my_axis = [-.9 .9 -.6 1.4];
inv_pc1 = 1;
inv_pc2 = 0;
figure('visible','off')
get_pca_mean_activity_mot_proj
end


if 0
get_legend;
end
