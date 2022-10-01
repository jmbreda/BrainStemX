%function Phases_polar_plot()
clear all; close all; clc;

set(0,'DefaultTextFontSize',24)
marker_size = 60;
% Saturation = r
% Hue = phi/2*pi

data_folder = '/scicore/home/zavolan/breda/BrainStemX/data/Bulk_Tables';
Tm = readtable([data_folder '/Gene_SampleName_MeanLogTPM_Pseudocount_CST_gt0.txt'],'ReadRowNames',1);
ind = find(~cellfun(@isempty,regexp(Tm.Properties.VariableNames,'NSC_WT')));
Tm = Tm(:,ind);
T = readtable([data_folder '/Gene_SampleName_LogTPM_Pseudocount_CST_gt0.txt'],'ReadRowNames',1);
ind = find(~cellfun(@isempty,regexp(T.Properties.VariableNames,'NSC_WT')));
T = T(:,ind);
tmp = cellfun(@(x) strsplit(x,'|'),T.Properties.RowNames,'UniformOutput',0);
tmp = [tmp{:}];
genes = tmp(2:2:end);



Am = Tm{:,:} - repmat(mean(Tm{:,:},2),1,size(Tm,2));
Am = Am';
A = T{:,:} - repmat(mean(T{:,:},2),1,size(T,2));
A = A';

samples = T.Properties.VariableNames; 
for s = 1:length(samples)
    tmp = strsplit(samples{s},'_');
    samples{s} = [tmp{1} '_' tmp{2} '_' tmp{3}];
end
C = load('../color_gradient_10.txt');
C = C/255;
my_grey = [.8 .8 .8];
Phases = {'Expansion','Neurogenesis','Gliogenesis'};

[~,~,idx_samp] = unique(samples,'stable');
my_color = C(idx_samp,:);

[coeff,score,latent,tsquared] = pca(Am);
% inverse pc 1 and 2
for pc = 1:2
    coeff(:,pc) = -coeff(:,pc);
    score(:,pc) = -score(:,pc);
end
% Get Projection of all samples in mean samples pc
score = A*coeff;

for p = Phases
    switch p{:}
    case 'Expansion'
        samp_i = 1:2;
    case 'Neurogenesis'
        samp_i = 3:7;
    case 'Gliogenesis'
        samp_i = 8:10;
    end
    idx.(p{:}) = zeros(size(idx_samp));
    for i = samp_i
        idx.(p{:}) = idx.(p{:}) + (idx_samp==i);
    end

    Score.(p{:}) = mean(score(find(idx.(p{:})),1:2),1);
    Score.(p{:}) = Score.(p{:})/norm(Score.(p{:}));
end

figure
hold on
vec_size = 125;
txt_size = 135;
for p = Phases
    plot(vec_size*[0 Score.(p{:})(1)],vec_size*[0 Score.(p{:})(2)],'-','color',my_grey);
    h = text(txt_size*Score.(p{:})(1),txt_size*Score.(p{:})(2),p{:});
    if Score.(p{:})(1)<0
        set(h,'HorizontalAlignment','right');
    end
end
scatter(score(:,1),score(:,2),marker_size,my_color,'o','filled')
xlabel(num2str(latent(1)/sum(latent),2))
ylabel(num2str(latent(2)/sum(latent),2))

% put Phi_0 (red) at Expansion
phi_0 =  atan2( Score.Expansion(1), Score.Expansion(2) );

R = linspace(125,130,8);
PHI = linspace(0,2*pi,401);
PHI = PHI(1:end-1);

x = [];
y = [];
col = [];
S = [];
for r = R
    for phi = PHI
        x = [x; r*sin(phi)];
        y = [y; r*cos(phi)];
        col = [col; hsv2rgb([ mod(phi-phi_0,2*pi)/(2*pi) r/max(R) 1])];
        S = [S; 5*r/max(R)];
    end
end
scatter(x,y,S.^2,col,'filled')

axis equal
set(gca,'XtickLabel',[],'YtickLabel',[])
axis off

if 1
    dim = [36 32];
    set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
    print(gcf,'Fig/Phase_plot_gene','-dpdf');
    print(gcf,'Fig/Phase_plot_gene','-dpng');
    print(gcf,'Fig/Phase_plot_gene','-dsvg');
end

% Get Phi for each gene
if 0
    phase_phi = zeros(3,1);
    for p = 1:length(Phases);
        phi_rad = atan2( Score.(Phases{p})(1), Score.(Phases{p})(2) );
        phase_phi(p,1) = mod(phi_rad-phi_0,2*pi)/(2*pi);
    end
    dlmwrite('phase_phi_gene.txt',phase_phi);

    phi = cell(length(genes),1);
    for m = 1:length(genes);
        phi_rad = atan2( coeff(m,1), coeff(m,2) );
        phi{m} = mod(phi_rad-phi_0,2*pi)/(2*pi);

        % Get phase
        dist = [abs(phase_phi-phi{m}), 1-abs(phase_phi-phi{m})];
        dist = min(dist')';
        [~,phase(m)] = min(dist);
    end
    gene_phi = containers.Map(genes,phi);
    gene_phase = containers.Map(genes,phase);
    clear gene;
    for p = 1:3
        my_genes{p} = genes(phase==p);
    end
    phase_gene = containers.Map([1:3],my_genes);

    save('gene_phi_NSC.mat','gene_phi');
    save('gene_phase_NSC.mat','gene_phase');
    save('phase_gene_NSC.mat','phase_gene');
end


