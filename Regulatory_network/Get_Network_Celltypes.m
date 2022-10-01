clear all; close all; clc;

celltype = 'NEUROGENESISuBPuNBN';
MY_MOT = {'Genes','All','Dlx1','E2f1','E2f2_E2f5','E2f8','Hey1_Myc_Mxi1','Hoxb7','Mef2c','Pax6','Sox3_Sox10','Tbx1_Eomes','Vsx1_Uncx_Prrx2_Shox2_Noto','Vsx2_Dlx3','Ybx1_Nfya_Nfyb_Nfyc_Cebpz','Zic3'};

MY_MOT = {'All'}
% Gene expression;
if 1
load('../../data/Matlab_Tables/Gene_SampleName_MeanLogTPM_Pseudocount_CST_gt2.mat');
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

for my_mot = MY_MOT(1)

    switch my_mot{:}
    case 'Genes'
        th_LL = 55;
    case 'All'
        th_LL = 30;
    case 'Ybx1_Nfya_Nfyb_Nfyc_Cebpz'
        th_LL = 50;
    otherwise
        th_LL = 20;
    end


    % Get Data :
    ismara_folder = ['/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna/' celltype '_WT'];
    [~,motifs] = textread([ismara_folder '/mat_indices'],'%u\t%s\n');
    sample =     textread([ismara_folder '/averaged_samples'],'%s\n')';
    genes =      textread([ismara_folder '/gene_indices'],'%s\n');
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


    if strcmp(my_mot{:},'All')
        % Get only genes with motif
        [genes,ind_gene,~] = intersect(genes,gene2motif.keys,'stable');
        Adj_m_g = Adj_m_g(:,ind_gene);

        % Get motifs with LL (in or out) higher than LL
        ind_mot = find(max(Adj_m_g')>th_LL);
        ind_gene = find(max(Adj_m_g)>th_LL);
    elseif strcmp(my_mot{:},'Genes')
        ind_mot = find(max(Adj_m_g')>th_LL);
        ind_gene = find(max(Adj_m_g)>th_LL);
    else
        % Get genes targeted by my motif and motif targeting my_motif
        [~,ind_my_mot,~] = intersect(motifs,my_mot,'stable');
        [~,ind_my_gene,~] = intersect(genes,strsplit(my_mot{:},'_'));

        ind_mot  = union(ind_my_mot, find(max(Adj_m_g(:,ind_my_gene),[],2)>th_LL));
        ind_gene = union(ind_my_gene,find(max(Adj_m_g(ind_my_mot,:),[],1)>th_LL));

        %Also add the targets of the motifs
        ind_gene = union(ind_gene,find(max(Adj_m_g(ind_mot,:),[],1)>th_LL));

    end
    
    % Add targeted genes as motifs
    if ~strcmp(my_mot{:},'Gene')
        for g = intersect(genes(ind_gene)',gene2motif.keys);
            ind_mot = union(ind_mot,find(strcmp(motifs,gene2motif(g{:}))));
        end
    end

    motifs = motifs(ind_mot);
    genes = genes(ind_gene);
    Adj_m_g = Adj_m_g(ind_mot,ind_gene);
    %Adj_m_m = Adj_m_m(ind_mot,:); 
    Activity = Activity(:,ind_mot);
    
    [tmp_genes,~,ind_E_gene] = intersect(genes,MeanLogBulk.Properties.RowNames,'stable');
    if length(tmp_genes) ~= length(genes)
        disp('!!! Genes not found in MeanLogBulk !!!')
    end
    ind_E_samp = intersect(sample,MeanLogBulk.Properties.VariableNames,'stable');
    GeneExp = MeanLogBulk(ind_E_gene,ind_E_samp);
    
    % Get A per phase normalized
    phase_ind.Expan = find(~cellfun(@isempty,regexp(sample,'NSC')));
    phase_ind.Neuro = find(~cellfun(@isempty,regexp(sample,'BP')));
    phase_ind.Glio  = find(~cellfun(@isempty,regexp(sample,'NBN')));

    k = 1;
    A = zeros(size(Activity,2),3);
    E = zeros(size(GeneExp,1),3);
    for p = {'Expan','Neuro','Glio'}
        A(:,k) = mean(Activity(phase_ind.(p{:}),:),1)';
        E(:,k) = mean(GeneExp{:,phase_ind.(p{:})},2);
        k = k+1;
    end
    A = (A - repmat(min(A')',1,3));
    A = A./repmat(sum(A,2),1,3);
    E = (E - repmat(min(E')',1,3));
    E = E./repmat(sum(E,2),1,3);

 
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
                    node.name = [node.name; genes{g}];
                    edge.target = [edge.target; genes{g}];
                    edge.gene = [edge.gene; ' '];
                else
                    node.name = [node.name; gene2motif(genes{g})];
                    edge.target = [edge.target; gene2motif(genes{g})]; 
                    edge.gene = [edge.gene; genes{g}];
                end
            end
        end
    end
    node.name = unique(node.name);

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
        if isempty(intersect(node.name{i},gene2motif.values))
            node.shape{i} = 'ellipse';
            node.image{i} = '';
            node.TF{i} = containers.Map();
            node.go{i} = {''};

            ind_gene = find(strcmp(node.name{i},genes));
            Phases_level = E(ind_gene,:);            
        else
            node.shape{i} = 'box';
            node.image{i} = [ismara_folder '/averaged_report/profiles/' node.name{i} '/profile_small.png'];
            corr_file = ['/scicore/home/zavolan/breda/BrainStemX/RegulatoryInteractions/correlation/' celltype '_WT/' node.name{i}];
            if exist(corr_file)
                [my_gene,my_corr] = textread(corr_file,'%s\t%f\n');
                node.TF{i} = containers.Map(my_gene,my_corr);
            else
                node.TF{i} = containers.Map();
            end

            ind_mot = find(strcmp(node.name{i},motifs));
            Phases_level = A(ind_mot,:);
            
            if isempty(intersect(node.name{i},motif2go.keys))
                node.go{i} = {' '};
            else
                node.go{i} = strsplit(motif2go(node.name{i}),';');
            end 
        end
        
        % get hsv color
        switch find(Phases_level(:,:)==0)
        case 1
            phi = Phases_level(:,2)*1/3 + Phases_level(:,3)*2/3;
        case 2
            phi = Phases_level(:,1)*1 + Phases_level(:,3)*2/3;
        case 3
            phi = Phases_level(:,1)*0 + Phases_level(:,2)*1/3;
        end
        hsv_color = [phi 1 1];
        rgb_color = hsv2rgb(hsv_color);
        node.color{i} = ['"' num2str(hsv_color(1)) ',' num2str(hsv_color(2)) ',' num2str(hsv_color(3)) '"'];

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
    
    % Graph properties
    graph_label = '<<table border="0"><tr><td><IMG SRC="Fig/Celltype_phase_plot.png"/></td></tr></table>>';
    if strcmp(my_mot,'All')
        rank_sep = '5';
    else
        rank_sep = '0';
    end
    %Graph = {'overlap=false','Goverlap=prism','fontname="Helvetica"',['label=' graph_label],['ranksep=' rank_sep],'ratio=compress'};
    Graph = {'overlap=false','Goverlap=prism','fontname="Helvetica"',['label=' graph_label,['ranksep=' rank_sep]],'rankdir=LR'};
    %Graph = {'overlap=false','Goverlap=prism','fontname="Helvetica"','rankdir=TB'};
    out_file = 'output/my_graph.dot';
    out_fig  = ['Fig/' celltype '/' my_mot{:}];
    
    % Plot Graph with graphviz
    Print_graphviz(Node,Edge,Graph,out_file,out_fig,my_mot)
end


