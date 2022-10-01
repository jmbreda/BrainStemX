%function Phases_polar_plot()
clear all; close all; clc;

set(0,'DefaultTextFontSize',24)
marker_size = 60;
% Saturation = r
% Hue = phi/2*pi

ismara_folder = '/scicore/home/zavolan/breda/BrainStemX/data/Ismara/Pseudocount_CST_gt2_no_mirna/NSC_WT';
[~,motifs] = textread([ismara_folder '/mat_indices'],'%s\t%s\n');
Am = load([ismara_folder '/averaged_activities']);
A = load([ismara_folder '/activities']);
samples = textread([ismara_folder '/sample_indices'],'%s');
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
vec_size = 1;
txt_size = 1.1;
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

R = linspace(1,1.05,8);
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
        S = [S; 5*r];
    end
end
scatter(x,y,S.^2,col,'filled')

axis equal
set(gca,'XtickLabel',[],'YtickLabel',[])
axis off

if 0
    dim = [36 32];
    set(gcf,'units','Centimeters','PaperUnits','Centimeters','PaperPositionMode','Auto','PaperPosition',[0 0 dim],'PaperSize',[dim]);
    print(gcf,'Fig/Phase_plot_motif','-dpdf');
    print(gcf,'Fig/Phase_plot_motif','-dpng');
end

% Get Phi for each motif
if 1
    phase_phi = zeros(3,1);
    for p = 1:length(Phases);
        phi_rad = atan2( Score.(Phases{p})(1), Score.(Phases{p})(2) );
        phase_phi(p,1) = mod(phi_rad-phi_0,2*pi)/(2*pi);
    end
    dlmwrite('phase_phi_motif.txt',phase_phi);

    phi = cell(length(motifs),1);
    for m = 1:length(motifs);
        phi_rad = atan2( coeff(m,1), coeff(m,2) );
        phi{m} = mod(phi_rad-phi_0,2*pi)/(2*pi);

        % Get phase
        dist = [abs(phase_phi-phi{m}), 1-abs(phase_phi-phi{m})];
        dist = min(dist')';
        [~,phase(m)] = min(dist);
    end
    motif_phi = containers.Map(motifs,phi);
    motif_phase = containers.Map(motifs,phase);
    clear my_motifs;
    for p=1:3
        my_motifs{p} = motifs(phase==p);
    end
    phase_motif = containers.Map([1:3],my_motifs);

    save('motif_phi_NSC.mat','motif_phi');
    save('motif_phase_NSC.mat','motif_phase');
    save('phase_motif_NSC.mat','phase_motif');
end


