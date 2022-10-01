clear all; close all; clc;


celltype = 'NSC';
MY_MOT = {'All'};
%MY_MOT = {'Chd1_Pml','Pou5f1','Hoxb7','Rest'};
%MY_MOT = {'Nfia','Nfic_Nfib','Nfix'};
%MY_MOT = {'Chd1_Pml','Pou5f1','Hoxb7','Rest','Nfia','Nfic_Nfib','Nfix'};
%MY_MOT = {'E2f1','Max_Mycn','Tead3_Tead4','Foxi1_Foxo1'};
%MY_MOT = {'Genes'};
%MY_MOT = {'Tead3_Tead4','Tead1'};

my_mot = strjoin(MY_MOT,'_');

switch my_mot
case 'All'
    th_LL = 28;
case 'Chd1_Pml_Pou5f1_Hoxb7_Rest'
    th_LL = 40;
case 'Nfia_Nfic_Nfib_Nfix'
    th_LL = 12;
case 'E2f1_Max_Mycn_Tead3_Tead4_Foxi1_Foxo1'
    th_LL = 20;
case 'Tead3_Tead4_Tead1'
    th_LL = 7;
otherwise
    th_LL = 30;
end

ismara_folder = ['/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna/' celltype '_WT'];

% Get Data :
[~,motifs] = textread([ismara_folder '/mat_indices'],'%u\t%s\n');
sample =     textread([ismara_folder '/averaged_samples'],'%s\n')';
genes =      textread([ismara_folder '/gene_indices'],'%s\n');
Adj_m_g =    load([ismara_folder '/Adjacency_motif_gene']);
%Adj_m_m =    load([ismara_folder '/Adjacency_motif_motif']);
load('gene2motif.mat');

load('motif_NSC_phi.mat');
phase_phi_motif = load('phase_phi_motif.txt');

load('gene_NSC_phi.mat');
phase_phi_gene = load('phase_phi_gene.txt');

!cut -f1 NSC_motifs_GO.txt > tmp_mot
!cut -f2 NSC_motifs_GO.txt > tmp_go
go_mot = textread('tmp_mot','%s','delimiter','\n');
go = textread('tmp_go','%s','delimiter','\n');
motif2go = containers.Map(go_mot,go);
!rm tmp*



if strcmp(my_mot,'All')
    % Get only genes with motif
    [genes,ind_gene,~] = intersect(genes,gene2motif.keys,'stable');
    Adj_m_g = Adj_m_g(:,ind_gene);

    % Get motifs with LL (in or out) higher than LL
    ind_mot = find(max(Adj_m_g')>th_LL);
    ind_gene = find(max(Adj_m_g)>th_LL);
elseif strcmp(my_mot,'Genes')
    [tot_LL,ind] = sort(sum(Adj_m_g,2),'descend');
    Adj_m_g = Adj_m_g(ind(1:40),:);
    motifs = motifs(ind(1:40));
    ind_mot = find(max(Adj_m_g')>th_LL);
    ind_gene = find(max(Adj_m_g)>th_LL);
else
    % Get genes targeted by my motif and motif targeting my_motif
    [~,ind_my_mot,~] = intersect(motifs,MY_MOT,'stable');
    [~,ind_my_gene,~] = intersect(genes,strsplit(my_mot,'_'));

    ind_mot  = union(ind_my_mot, find(max(Adj_m_g(:,ind_my_gene),[],2)>th_LL));
    ind_gene = union(ind_my_gene,find(max(Adj_m_g(ind_my_mot,:),[],1)>th_LL));

    [~,ind_gene_mot,~] = intersect(genes,gene2motif.keys,'stable');
    ind_gene = intersect(ind_gene,ind_gene_mot);
end

% Add targeted genes as motifs
if ~strcmp(my_mot,'Gene')
    for g = intersect(genes(ind_gene)',gene2motif.keys);
        ind_mot = union(ind_mot,find(strcmp(motifs,gene2motif(g{:}))));
    end
end

motifs = motifs(ind_mot);
genes = genes(ind_gene);
Adj_m_g = Adj_m_g(ind_mot,ind_gene);
%Adj_m_m = Adj_m_m(ind_mot,:); 


% Get A per phase normalized
phase_ind.Expan = find(~cellfun(@isempty,regexp(sample,'E1[0-1]5')));
phase_ind.Neuro = find(~cellfun(@isempty,regexp(sample,'E1[2-6]5')));
phase_ind.Glio  = [find(~cellfun(@isempty,regexp(sample,'E1[7-8]5'))) find(~cellfun(@isempty,regexp(sample,'PN')))];

% Initiallize node and edges
edge.source = {};
edge.target = {};
edge.ll = [];
edge.gene = {};
node.name = {};
% Get Motif -> Gene egdes : Loop over marker genes
%[M,G] = fidt
for m = 1:length(motifs)
    node.name = [node.name; motifs{m}];
    for g = 1:length(genes)
        if Adj_m_g(m,g)>th_LL
            edge.ll = [edge.ll; Adj_m_g(m,g)];
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

% Make Edge Table
Edge = table(edge.source,edge.target,edge.ll,edge.gene,'VariableNames',{'source','target','ll','gene'});

% Make Gene profile
%Plot_gene_profile(intersect(node.gene,Gene_log_tpm.Gene));

% Initiate node attributes
node.shape = cell(length(node.name),1);
node.image = cell(length(node.name),1);
node.color = cell(length(node.name),1);
node.phase = zeros(length(node.name),1);
node.go = cell(length(node.name),1);
node.TF = cell(length(node.name),1);
node.in_degree = zeros(length(node.name),1);
node.out_degree = zeros(length(node.name),1);


% Make motif nodes
for i = 1:length(node.name)
    if isempty(intersect(node.name{i},gene2motif.values))
        node.shape{i} = 'ellipse';
        node.image{i} = '';
        node.TF{i} = containers.Map();
        node.go{i} = {''};

        node.color{i} = ['"' num2str(gene_phi(node.name{i})) ',1,1"'];

        % Get phase of my motif
        dist = [abs(phase_phi_gene-gene_phi(node.name{i})), 1-abs(phase_phi_gene-gene_phi(node.name{i}))];
        dist = min(dist')';
        [~,node.phase(i)] = min(dist);
    else
        node.shape{i} = 'box';
        node.image{i} = '';%['profile/Motif/NSC_WT/png/' node.name{i} '_small.png'];
        corr_file = ['/scicore/home/zavolan/breda/BrainStemX/RegulatoryInteractions/correlation/NSC_WT/' node.name{i}];
        if exist(corr_file)
            [my_gene,my_corr] = textread(corr_file,'%s\t%f\n');
            node.TF{i} = containers.Map(my_gene,my_corr);
        else
            node.TF{i} = containers.Map();
        end

        if isempty(intersect(node.name{i},motif2go.keys))
            node.go{i} = {''};
        else
            node.go{i} = strsplit(motif2go(node.name{i}),';');
        end
    
        node.color{i} = ['"' num2str(motif_phi(node.name{i})) ',1,1"'];

        % Get phase of my motif
        dist = [abs(phase_phi_motif-motif_phi(node.name{i})), 1-abs(phase_phi_motif-motif_phi(node.name{i}))];
        dist = min(dist')';
        [~,node.phase(i)] = min(dist);
    end
    node.in_degree(i) = sum(strcmp(edge.target,node.name{i}));
    node.out_degree(i) = sum(strcmp(edge.source,node.name{i}));
    
end

% Make Node table
Node = table(node.name,node.shape,node.image,node.TF,node.color,node.go,node.phase,node.in_degree,node.out_degree,'VariableNames',{'name','shape','image','TF','color','go','phase','in_degree','out_degree'});

% Graph properties
%graph_label = '<<table border="0"><tr><td><IMG SRC="Fig/Phase_plot_motif.png"/></td></tr></table>>';
if strcmp(my_mot,'All')
    rank_sep = '20';
else
    rank_sep = '0';
end
%Graph = {'overlap=false','Goverlap=prism','fontname="Helvetica"',['label=' graph_label],['ranksep=' rank_sep],'rankdir=LR'};
Graph = {'overlap=false','Goverlap=prism','fontname="Helvetica"',['ranksep=' rank_sep],'rankdir=LR','concentrate=false'};
out_file = 'output/my_graph.dot';
out_fig  = ['Fig/' celltype '/' my_mot];

% Plot Graph with graphviz
Print_graphviz_NSC(Node,Edge,Graph,out_file,out_fig,my_mot)


