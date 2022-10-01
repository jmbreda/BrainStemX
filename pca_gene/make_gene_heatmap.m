clc; close all; clear all;

load('my_data')

subplot_line = 19;
subplot_col  = 9;
subplot_index = reshape(1:(subplot_line*subplot_col),subplot_col,subplot_line)';


default_fs = 6;
set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',default_fs)
dim = [17 26];
my_grey = [.8 .8 .8];

%figure('visible','off')


% PCA ALL orth to NSC
my_title = 'Celltypes';
l_proj = 900;
n_best = 16;
n_clust = 3;
my_axis = [-55 40 -40 60];
inv_pc1 = 0;

Line = 5;
Col = 0;
%get_subplot_gene_pca_All


% PCA NSC
celltype = 'NSC'
my_title = 'Neural stem cells';
my_marker = 'o';
n_genes = 50;

Line = 5;
Col = 0;
get_gene_heatmap_NSC_BP_NBN


% Print figure
if 0
set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0   dim],'PaperSize',[dim]);
print(gcf,'Fig/Gene_dynamic','-dpdf');
end
