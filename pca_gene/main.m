clc; close all; clear all;

addpath('~/Matlab_Scripts/cbrewer/')

load('my_data')

default_fs = 6;
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',default_fs)
dim = [8 6];
my_grey = [.8 .8 .8];

Timepoint = {'E105','E115','E125','E135','E145','E155','E165','E175','E185','PN'};
%C = cbrewer('div','Spectral',10);
C = cbrewer('div','RdYlGn',10);

CellTypes = {'NSC','BP','NBN'};
Markers = {'s','o','^'};


if 1
% PCA ALL orth to NSC
my_title = 'Celltypes';
celltype = 'ALL'
l_proj = 900;
n_best = 20;
n_clust = 3;
my_axis = [-50 50 -40 100];
inv_pc1 = 0;
inv_pc2 = 0;
get_subplot_gene_pca_All
end

n_best = 20;

if 1
celltype = 'NEUROGENESISuBPuNBN'
l_proj = 1600;
n_clust = 6;
my_axis = [-90 80 -130 120];
inv_pc1 = 0;
inv_pc2 = 0;
get_subplot_gene_pca_NSC_BP_NBN
end

if 1
% PCA NSC
celltype = 'NSC'
l_proj = 1150;
n_clust = 5;
my_axis = [-120 80 -60 120];
inv_pc1 = 1;
inv_pc2 = 0;
get_subplot_gene_pca_NSC_BP_NBN
end

if 1
% PCA BP
celltype = 'BP'
l_proj = 800;
n_clust = 5;
my_axis = [-65 60 -60 60];
inv_pc1 = 0;
inv_pc2 = 0;
get_subplot_gene_pca_NSC_BP_NBN
end

if 1
% PCA NBN
celltype = 'NBN'
l_proj = 500;
n_clust = 5;
my_axis = [-40 40 -100 40];
inv_pc1 = 1;
inv_pc2 = 1;
get_subplot_gene_pca_NSC_BP_NBN
end
