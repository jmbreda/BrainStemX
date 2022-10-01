clear all; close all; clc;


celltype = 'NSC';
Phases = {'Expansion','Neurogenesis','Gliogenesis'};

for p = 1:3
    switch p
    case 1
        th_LL = 30;
    case 2
        th_LL = 35;
    case 3
        th_LL = 25;
    end

    ismara_folder = ['/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna/' celltype '_WT'];

    % Get Data :
    [~,motifs] = textread([ismara_folder '/mat_indices'],'%u\t%s\n');
    sample =     textread([ismara_folder '/averaged_samples'],'%s\n')';
    genes =      textread([ismara_folder '/gene_indices'],'%s\n');
    Adj_m_g =    load([ismara_folder '/Adjacency_motif_gene']);
    load('gene2motif.mat');

    load('motif_phi_NSC.mat');
    load('motif_phase_NSC.mat');
    load('phase_motif_NSC.mat');

    load('gene_phi_NSC.mat');
    load('gene_phase_NSC.mat');
    load('phase_gene_NSC.mat');

    !cut -f1 motifs_GO.txt > tmp_mot
    !cut -f2 motifs_GO.txt > tmp_go
    go_mot = textread('tmp_mot','%s','delimiter','\n');
    go = textread('tmp_go','%s','delimiter','\n');
    motif2go = containers.Map(go_mot,go);
    !rm tmp*

    [~,ind_mot,~] = intersect(motifs,phase_motif(p));
    [~,ind_gene,~] = intersect(genes,phase_gene(p));
    motifs = motifs(ind_mot);
    genes = genes(ind_gene);
    Adj_m_g = Adj_m_g(ind_mot,ind_gene);

    ind_mot = find(max(Adj_m_g')>th_LL);
    ind_gene = find(max(Adj_m_g)>th_LL);
    % Add targeted genes as motifs
    for g = intersect(genes(ind_gene)',gene2motif.keys);
        ind_mot = union(ind_mot,find(strcmp(motifs,gene2motif(g{:}))));
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
            node.phase(i) = gene_phase(node.name{i});
        else
            node.shape{i} = 'box';
            node.image{i} = ['profile/Motif/' celltype '_WT/' node.name{i} '_small.png'];
            corr_file = ['/scicore/home/zavolan/breda/BrainStemX/RegulatoryInteractions/correlation/' celltype '_WT/' node.name{i}];
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
            node.phase(i) = motif_phase(node.name{i}); 
        end
        node.in_degree(i) = sum(strcmp(edge.target,node.name{i}));
        node.out_degree(i) = sum(strcmp(edge.source,node.name{i}));
    end

    % Make Node table
    Node = table(node.name,node.shape,node.image,node.TF,node.color,node.go,node.phase,node.in_degree,node.out_degree,'VariableNames',{'name','shape','image','TF','color','go','phase','in_degree','out_degree'});

    % Graph properties
    rank_sep = '8';
    Graph = {'overlap=false','Goverlap=prism','fontname="Helvetica"',['ranksep=' rank_sep],'rankdir=LR','concentrate=false'};
    out_file = 'output/my_graph.dot';
    out_fig  = ['Fig/' celltype '_' Phases{p}];

    % Plot Graph with graphviz
    Print_graphviz_NSC_Phases(Node,Edge,Graph,out_file,out_fig)
end

