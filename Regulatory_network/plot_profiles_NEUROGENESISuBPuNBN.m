function plot_profiles_NEUROGENESISuBPuNBN()

    purple = [0.4940 0.1840 0.5560];
    green = [0.4660 0.6740 0.1880];
    blue = [0.3010 0.7450 0.9330];

    celltype = 'NEUROGENESISuBPuNBN_WT'

    for data = {'Gene'}% 'Motif'}
        disp(data{:})
    
        %mkdir(['profile/' data{:} '/' celltype])

        if strcmp(data{:},'Gene');
            load('~/BrainStemX/data/Matlab_Tables/Gene_SampleName_MeanLogTPM_Pseudocount_CST_gt0.mat')
            load('~/BrainStemX/data/Matlab_Tables/Gene_SampleName_StdDevLogTPM_Pseudocount_CST_gt0.mat')

            idx_samp = [3:7 11:23]+1;
            Samples = MeanLogBulk.Properties.VariableNames(idx_samp);

            my_genes = {'Neurog1','Hmga2','Sox3','Hbb-y','Crabp2','Fezf1','Dhrs3','Nfyb','Pou3f1','E2f1','Tbp','Foxo1','Foxd1'};

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
            plot_profile_small(E(i,:),['profile/' data{:} '/' celltype '/' Name{i} '_small'])
        end
    end
end

function plot_profile(E,dE,name,file,data)
    default_fs = 24;
    set(0,'DefaultAxesFontName','Arial','DefaultAxesFontSize',default_fs)
    
    idx.nsc = 1:5;
    color.nsc = [0.1059    0.6196    0.4667];
    t.nsc = 12.5:16.5;
    
    idx.bp = 6:13;
    color.bp = [0.8510    0.3725    0.0078];
    t.bp = 12.5:19.5;
    
    idx.nbn = 14:18;
    color.nbn = [0.4588    0.4392    0.7020];
    t.nbn = 15.5:19.5;

    figure('visible','off')
    hold on 

    for ct = {'nsc','bp','nbn'}
        errorbar(t.(ct{:}),E(idx.(ct{:})),dE(idx.(ct{:})),'color',color.(ct{:}),'linewidth',2) 
    end
    title(name,'FontSize',32,'interpreter','none')
    if strcmp(data,'Motif')
        ylabel('activity')
    else
        ylabel('log expression')
    end
    xlabel('time')
    axis([12 20 min(E-2*dE) max(E+2*dE)]);
    set(gca,'Xtick',[],'YTick',[])
    dim_pxl = [400 300];
    dim_pnt = dim_pxl*.75;
    %set(gca,'position',[.1 . .9 .9]);
    set(gcf,'units','pixels','Position',[0 0 dim_pxl],...
            'PaperUnits','points','PaperPositionMode','Auto',...
            'PaperPosition',[0 0 dim_pnt],'PaperSize',[dim_pnt])
    print(gcf,file,'-dsvg');
end

function plot_profile_small(E,file)
    
	idx.nsc = 1:5;
    color.nsc = [0.1059    0.6196    0.4667];
    t.nsc = 12.5:16.5;
    
    idx.bp = 6:13;
    color.bp = [0.8510    0.3725    0.0078];
    t.bp = 12.5:19.5;
    
    idx.nbn = 14:18;
    color.nbn = [0.4588    0.4392    0.7020];
    t.nbn = 15.5:19.5;

   
    figure('visible','off')
    hold on
	for ct = {'nsc','bp','nbn'}
        plot(t.(ct{:}),E(idx.(ct{:})),'o-','color',color.(ct{:}),'MarkerEdgeColor',[1 1 .99],'MarkerFaceColor',color.(ct{:}),'linewidth',2,'MarkerSize',8)
    end
    axis([12 20 1.05*min(E) 1.05*max(E)]);
    axis off;
    dim_pxl = [400 200];
    dim_pnt = dim_pxl*.75;
    set(gca,'position',[.05 .05 .9 .9]);
    set(gcf,'units','pixels','Position',[0 0 dim_pxl],...
            'PaperUnits','points','PaperPositionMode','Auto',...
            'PaperPosition',[0 0 dim_pnt],'PaperSize',[dim_pnt])
    print(gcf,file,'-dpng');
end


