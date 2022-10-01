function plot_profiles()

    celltypes = {'NSC_WT'};

    for data = {'Gene','Motif'}
        disp(data{:})
        for c = 1:length(celltypes)
            celltype = celltypes{c};
            disp(celltype)

            mkdir(['profile/' data{:} '/' celltype])

            if strcmp(data{:},'Gene');
                load('~/BrainStemX/data/Matlab_Tables/Gene_SampleName_MeanLogTPM_Pseudocount_CST_gt0.mat')
                load('~/BrainStemX/data/Matlab_Tables/Gene_SampleName_StdDevLogTPM_Pseudocount_CST_gt0.mat')

                idx_samp = find(~cellfun(@isempty,regexp(MeanLogBulk.Properties.VariableNames,'NSC_WT')));
                Samples = MeanLogBulk.Properties.VariableNames(idx_samp);

                my_genes = {'Neurod2','Neurod6','Nfix','Crabp2','Hmga2','Hbb-y','Hba-x','Hbb-bh1','Tnc','Slco1c1','Bcan','Scrg1','Aldh1l1','Luzp2','Olig1','Hepacam','Sparcl1','Gfap','Slc6a11','Pdgfra','Hoxb7','Neurod1','Tef','Pou5f1','Sox5','Sox10','Sox2','Fabp7','Nfia','Nfic','Stat2','Tgif1','E2f1','Max'};

                E = zeros(length(my_genes),length(idx_samp));
                dE = E;
                for g = 1:length(my_genes)
                    idx_gene = find(~cellfun(@isempty,regexp(MeanLogBulk.GeneID,my_genes{g})));
                    E(g,:) = sum(MeanLogBulk{idx_gene,idx_samp},1);
                    dE(g,:) = mean(StdLogBulk{idx_gene,idx_samp},1);
                end
                Name = my_genes;

            elseif strcmp(data{:},'Motif')

                ismara_folder = ['/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna/' celltype];
                [~,Name] = textread([ismara_folder '/mat_indices'],'%u\t%s\n');
                E = load([ismara_folder '/averaged_activities'])';
                dE = load([ismara_folder '/averaged_deltas'])';

            end

            close all
            parfor i = 1:length(Name)'
                plot_profile(E(i,:),dE(i,:),Name{i},['profile/' data{:} '/' celltype '/' Name{i}],data{:})
                %plot_profile_small(E(i,:),purple,['profile/' data{:} '/' celltype '/' Name{i} '_small'])
            end
        end
    end
end

function plot_profile(E,dE,name,file,data)
    default_fs = 24;
    set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',default_fs)
    
    color = [0.1059    0.6196    0.4667];

    figure('visible','off')
    hold on 

    t = 1:length(E);
    errorbar(t,E,dE,'color',color,'linewidth',2)
    title(name,'FontSize',32,'interpreter','none')
    if strcmp(data,'Motif')
        ylabel('activity')
    else
        ylabel('log expression')
    end
    xlabel('time')
    axis([min(t)-1 max(t)+1 min(E-2*dE) max(E+2*dE)]);
    set(gca,'Xtick',[],'YTick',[])
    dim_pxl = [400 300];
    dim_pnt = dim_pxl*.75;
    %set(gca,'position',[.1 . .9 .9]);
    set(gcf,'units','pixels','Position',[0 0 dim_pxl],...
            'PaperUnits','points','PaperPositionMode','Auto',...
            'PaperPosition',[0 0 dim_pnt],'PaperSize',[dim_pnt])
    print(gcf,file,'-dsvg');
end

function plot_profile_small(E,color,file)
    
    figure('visible','off')
    hold on 

    t = 1:length(E);
    plot(t,E,'o-','color',color,'MarkerEdgeColor',[1 1 .99],'MarkerFaceColor',color,'linewidth',2,'MarkerSize',8)
    axis([min(t) max(t) 1.05*min(E) 1.05*max(E)]);
    axis off;
    dim_pxl = [200 100];
    dim_pnt = dim_pxl*.75;
    set(gca,'position',[.05 .05 .9 .9]);
    set(gcf,'units','pixels','Position',[0 0 dim_pxl],...
            'PaperUnits','points','PaperPositionMode','Auto',...
            'PaperPosition',[0 0 dim_pnt],'PaperSize',[dim_pnt])
    print(gcf,file,'-dpng');
end
