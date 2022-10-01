clear all; close all; clc;

addpath('~/Matlab_Scripts/cbrewer/');

celltype2time = cbrewer('div','Spectral',128);
celltype2time = cool(128);
% Get legend
if 1
	figure('visible','off');
	imagesc(linspace(0,1,128));
	colormap(celltype2time);
	axis off;
	dim = [20 4];
	set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0  0 dim],'PaperSize',[dim]);
	print(gcf,'Fig/NEUROGENESISuBPuNBN/color_legend','-dpng');	
	print(gcf,'Fig/NEUROGENESISuBPuNBN/color_legend','-dsvg');
end

% Motif activity :
ismara_folder = ['/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna/NEUROGENESISuBPuNBN_WT'];
[~,Motifs] = textread([ismara_folder '/mat_indices'],'%u\t%s\n');
Sample =     textread([ismara_folder '/averaged_samples'],'%s\n')';
Genes =      textread([ismara_folder '/gene_indices'],'%s\n');
Activity =   load([ismara_folder '/averaged_activities']);
Adj_m_g =    load([ismara_folder '/Adjacency_motif_gene']);
%Adj_m_m =    load([ismara_folder '/Adjacency_motif_motif']);
load('gene2motif.mat');

!cut -f1 NEUROGENESISuBPuNBN_motifs_GO.txt > tmp_mot
!cut -f2 NEUROGENESISuBPuNBN_motifs_GO.txt > tmp_go
go_mot = textread('tmp_mot','%s','delimiter','\n');
go = textread('tmp_go','%s','delimiter','\n');
motif2go = containers.Map(go_mot,go);
!rm tmp*

[keys,vals] = textread([ismara_folder '/active_matrices'],'%s\t%f\n');
motif2Zval = containers.Map(keys,vals);
clear keys vals;

[coeff,score,latent] = pca(Activity);


% Gene expression;
if 1 
    load('../../data/Matlab_Tables/Gene_SampleName_MeanLogTPM_Pseudocount_CST_gt2.mat');

    % Get NEUROGENESISuBPuNBN WT samples
    idx_ko = find(~cellfun(@isempty,regexp(MeanLogBulk.Properties.VariableNames,'KO')));
    idx_Expan_Glio = find(~cellfun(@isempty,regexp(MeanLogBulk.Properties.VariableNames,'E105_NSC|E115_NSC|E175_NSC|E185_NSC|PN_NSC')));
    MeanLogBulk(:,union(idx_ko,idx_Expan_Glio)) = [];

    tmp = cellfun(@(x) strsplit(x,'|'),MeanLogBulk.GeneID,'UniformOutput',0);
    tmp = [tmp{:}];
    genes = tmp(2:2:end);
    samp = MeanLogBulk.Properties.VariableNames(2:end);

    M = MeanLogBulk{:,2:end};
    [genes,~,uniq_ind] = unique(genes);
    [xx, yy] = ndgrid(uniq_ind,1:size(M,2));
    C=accumarray([xx(:) yy(:)],M(:));

    MeanLogBulk = array2table(C);
    MeanLogBulk.Properties.VariableNames = samp;
    MeanLogBulk.Properties.RowNames = genes;

    clear M C genes samp xx yy tmp;
end


% Get A per phase normalized
phase_ind.NSC = find(~cellfun(@isempty,regexp(Sample,'NSC')));
phase_ind.BP = find(~cellfun(@isempty,regexp(Sample,'BP')));
phase_ind.NBN  = find(~cellfun(@isempty,regexp(Sample,'NBN')));

My_Component = {'celltype','time','both'};
for cmp = 1:3
	
	if cmp==1 || cmp==2 || cmp==3
		switch cmp
		case 1
			th_LL = 15;
			th_pc = 1e-3;
		case 2
			th_LL = 15;
			th_pc = 1e-3;
		case 3
			th_LL = 15;
			th_pc = 0;
		end


		% Get only genes with motif
		[~,idx_gene,~] = intersect(Genes,gene2motif.keys,'stable');

		% Get only motifs above pca threshold
		if cmp==3
			idx_mot = find( sum(coeff(:,1:2).^2,2) > 0 );
		else
			idx_mot = find( coeff(:,cmp).^2 > th_pc );
		end

		%[~,idx_totvar]=sort(coeff(:,1:2).^2*latent(1:2),'descend');
		%idx_mot = intersect(idx_mot,idx_totvar(1:20));

		% Get motifs and genes with LL higher than th_LL
		idx_mot = intersect(idx_mot,find(sum(Adj_m_g')>th_LL));
		idx_gene = intersect(idx_gene,find(sum(Adj_m_g)>th_LL));
		
		% Add genes as potential motif
		idx_mot_from_gene = [];
		for i = idx_gene'
			idx_mot_from_gene(end+1) = find(strcmp(Motifs,gene2motif(Genes{i})));
		end

		idx_gene_from_mot = [];
		for i = idx_mot'
			for g = strsplit(Motifs{i},'_')
				idx_tmp = find(strcmp(Genes,g));
				if ~isempty(idx_tmp)
					idx_gene_from_mot(end+1) = idx_tmp;
				end
			end
		end

		%idx_mot = union(idx_mot,idx_mot_from_gene);
		%idx_gene = union(idx_gene,idx_gene_from_mot);

		% Get subsets
		motifs = Motifs(idx_mot);
		genes = Genes(idx_gene);
		adj_m_g = Adj_m_g(idx_mot,idx_gene);

		% Initiallize node and edges
		edge.source = {};
		edge.target = {};
		edge.ll = [];
		edge.gene = {};
		% Get Motif -> Gene egdes
		for m = 1:length(motifs)
			for g = 1:length(genes)
				if adj_m_g(m,g)>th_LL && ( motif2Zval(motifs{m}) > 1 ) && ( motif2Zval(gene2motif(genes{g})) > 1 )
					edge.ll = [edge.ll; adj_m_g(m,g)];
					edge.source = [edge.source; motifs{m}];

					if isempty(intersect(genes{g},gene2motif.keys))
						edge.target = [edge.target; genes{g}];
						edge.gene = [edge.gene; ' '];
					else
						edge.target = [edge.target; gene2motif(genes{g})]; 
						edge.gene = [edge.gene; genes{g}];
					end
				end
			end
		end
		node.name = union(edge.source,edge.target);

		% Make sure edge name doesn't start with a number
		for i = 1:length(edge.target)
			if ~isstrprop(edge.target{i}(1),'alpha')
				edge.target{i} = ['_' edge.target{i}];
			end
		end
		edge.source = strrep(edge.source,'-','_');
		edge.target = strrep(edge.target,'-','_');
		edge.source = strrep(edge.source,'.','_');
		edge.target = strrep(edge.target,'.','_');

		% Make Edge Table
		Edge = table(edge.source,edge.target,edge.ll,edge.gene,'VariableNames',{'source','target','ll','gene'});
		Edge_both.(My_Component{cmp}) = Edge;

		% Make Gene profile
		%Plot_gene_profile(intersect(node.gene,Gene_log_tpm.Gene));

		% Initiate node attributes
		node.shape = cell(length(node.name),1);
		node.image = cell(length(node.name),1);
		node.color = cell(length(node.name),1);
		node.go = cell(length(node.name),1);
		node.TF = cell(length(node.name),1); 

		% Make motif nodes
		for i = 1:length(node.name)
			node.shape{i} = 'box';
			node.image{i} = '';%['profile/Motif/NEUROGENESISuBPuNBN_WT/' node.name{i} '_small.png'];
			corr_file = ['/scicore/home/zavolan/breda/BrainStemX/RegulatoryInteractions/correlation/NEUROGENESISuBPuNBN_WT/' node.name{i}];
			if exist(corr_file)
				[my_gene,my_corr] = textread(corr_file,'%s\t%f\n');
				node.TF{i} = containers.Map(my_gene,my_corr);
			else
				node.TF{i} = containers.Map();
			end
		
			if isempty(intersect(node.name{i},motif2go.keys))
				node.go{i} = {' '};
			else
				node.go{i} = strsplit(motif2go(node.name{i}),';');
			end 
		
			idx_mot_tmp = find(strcmp(node.name{i},Motifs));
			pc = coeff(idx_mot_tmp,1:2).^2;	   
			ct_score = pc(1)/sum(pc);
			my_color = rgb2hsv(celltype2time( round( 1 + ct_score*(length(celltype2time) - 1) ),:));
			node.color{i} = ['"' num2str(my_color(1)) ',' num2str(my_color(2)) ',' num2str(my_color(3)) '"'];

		end
		for i = 1:length(node.name)
			if ~isstrprop(node.name{i}(1),'alpha')
				node.name{i} = ['_' node.name{i}];
			end
		end
		node.name = strrep(node.name,'.','_');
		node.name = strrep(node.name,'-','_');

		% Make Node table
		Node = table(node.name,node.shape,node.image,node.TF,node.color,node.go,'VariableNames',{'name','shape','image','TF','color','go'});
		Node_both.(My_Component{cmp}) = Node;

	elseif cmp==3
		Edge = [Edge_both.time; Edge_both.celltype];
		Node = [Node_both.time; Node_both.celltype];
		[~,idx,~] = unique(Node.name);
		Node = Node(idx,:);
	end
    
    % Graph properties
    graph_label = '<<table border="0"><tr><td><IMG SRC="Fig/NEUROGENESISuBPuNBN/color_legend.png"/></td></tr></table>>';
    rank_sep = '0';
    %Graph = {'overlap=false','Goverlap=prism','fontname="Helvetica"',['label=' graph_label],['ranksep=' rank_sep],'ratio=compress'};
    Graph = {'overlap=false','Goverlap=prism','fontname="Helvetica"',['label=' graph_label,['ranksep=' rank_sep]],'rankdir=LR'};
    %Graph = {'overlap=false','Goverlap=prism','fontname="Helvetica"','rankdir=TB'};
    out_file = 'output/my_graph.dot';
    out_fig  = ['Fig/NEUROGENESISuBPuNBN/' My_Component{cmp}];
    
    % Plot Graph with graphviz
    Print_graphviz(Node,Edge,Graph,out_file,out_fig,'All')
end


