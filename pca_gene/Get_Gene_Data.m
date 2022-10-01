function Get_Gene_Data()

    Timepoint = {'E105','E115','E125','E135','E145','E155','E165','E175','E185','PN'};
    C = load('../color_gradient_10.txt');
    C = C/255;
    my_marker.NSC = 'o';
    my_marker.BP = '^';
    my_marker.NBN = 'h';

    load('~/BrainStemX/data/Matlab_Tables/Gene_SampleName_LogTPM_Pseudocount_CST_gt2.mat');
    load('~/BrainStemX/data/Matlab_Tables/Gene_SampleName_MeanLogTPM_Pseudocount_CST_gt2.mat');
    load('~/BrainStemX/data/Matlab_Tables/Gene_SampleName_StdDevLogTPM_Pseudocount_CST_gt2.mat');

    Tab.All = LogTPM{:,2:end};
    Tab.Mean = MeanLogBulk{:,2:end};
    Tab.Std = StdLogBulk{:,2:end};
    Genes = LogTPM.GeneID;
    g = cellfun(@(s) strsplit(s,'|'),Genes,'UniformOutput',0);
    g = [g{:}];
    Genes = g(2:2:end);
    Samples.All = LogTPM.Properties.VariableNames(2:end);
    Samples.Mean = MeanLogBulk.Properties.VariableNames(2:end);
    Samples.Std = StdLogBulk.Properties.VariableNames(2:end);

    for ii = fieldnames(Samples)'
        for ct = {'NSC','BP','NBN'}
            celltype_ind.(ii{:}).(ct{:}) = find( ~cellfun(@isempty,regexp(Samples.(ii{:}),[ct{:} '_WT'])) );
        end
        celltype_ind.(ii{:}).NEUROGENESISuBPuNBN = union(celltype_ind.(ii{:}).BP,celltype_ind.(ii{:}).NBN);
        NEUROG = find( ~cellfun(@isempty,regexp(Samples.(ii{:}),['E1[2-6]5_NSC_WT'])) );
        celltype_ind.(ii{:}).NEUROGENESISuBPuNBN = union(celltype_ind.(ii{:}).NEUROGENESISuBPuNBN,NEUROG);
    end

    save('my_data')
