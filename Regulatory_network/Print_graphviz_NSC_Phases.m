function Make_Graph_for_graphviz(Node,Edge,Graph,out_file,out_fig)

    disp('Plotting Graph')

    %brain_layer_marker = textread('brain_layer_marker.txt','%s','delimiter','\n');

    font_size = '32';
    linewidth = '14';
    arrow_size = '3';
    min_len = '1';

    N_Edge = size(Edge,1);

    for i = 1:length(Edge.target) 
        if ~isstrprop(Edge.target{i}(1),'alpha') 
            Edge.target{i} = ['_' Edge.target{i}]; 
        end 
    end
    for i = 1:length(Node.name) 
        if ~isstrprop(Node.name{i}(1),'alpha') 
            Node.name{i} = ['_' Node.name{i}];
        end 
    end
    Edge.source = strrep(Edge.source,'-','_');
    Edge.source = strrep(Edge.source,'.','_');
    Edge.target = strrep(Edge.target,'-','_');
    Edge.target = strrep(Edge.target,'.','_');
    Node.name = strrep(Node.name,'-','_');
    Node.name = strrep(Node.name,'.','_');


    fid = fopen(out_file,'w');

    fprintf(fid,'digraph G {\n');

    % Graph options
    fprintf(fid,'\tgraph[ ');
    for i = 1:length(Graph)
        fprintf(fid,'\t\t%s,\n',Graph{i});
    end
    fprintf(fid,'\t];\n');
    
    % Create Edges
    for i = 1:N_Edge
        arrow_head = 'normal';
        tar_ind = find(strcmp(Node.name,Edge.target{i}));
        if Node.TF{tar_ind}.Count > 0
            node_corr = Node.TF{tar_ind}(Edge.gene{i});
            if sign(node_corr)==-1
                arrow_head = 'inv';
            end
        end
        %fprintf(fid,'\t%s -> %s [arrowsize="%s",arrowhead="%s",fontsize=%s,minlen=%s,headlabel="%s",taillabel="%s"];\n',Edge.source{i},Edge.target{i},arrow_size,arrow_head,font_size,min_len,Edge.source{i},Edge.gene{i});
        fprintf(fid,'\t%s -> %s [arrowsize="%s",arrowhead="%s",fontsize=%s,minlen=%s];\n',Edge.source{i},Edge.target{i},arrow_size,arrow_head,font_size,min_len);
    end

    % Create Nodes
    for i = 1:size(Node,1);
        if strcmp(Node.shape{i},'ellipse')
            label = Node.name{i};
        else
            label = ['<<table border="0" cellspacing="0">'];
            label = [label '<tr><td>'];
            if Node.TF{i}.Count > 0
                for tf = Node.TF{i}.keys
                    label = [label tf{:} ':' num2str(Node.TF{i}(tf{:})) ' | '];
                end
                label(end-2:end) = '';
            else
                label = [label Node.name{i}];
            end
            label = [label '</td></tr>'];
            label = [label '<tr><td> </td></tr>'];
            label = [label '<tr><td> </td></tr>'];
            label = [label '<tr><td> </td></tr>'];
            for go = Node.go{i}
                label = [label '<tr><td>' go{:} '</td></tr>'];
            end
            label = [label '</table>>'];
        end

        fprintf(fid,'\t\t%s [shape=%s,label=%s,image="%s",color=%s,style="setlinewidth(%s)",fontsize=%s];\n',Node.name{i},Node.shape{i},label,Node.image{i},Node.color{i},linewidth,font_size);
    end
    fprintf(fid,'}\n');
    fclose(fid);

    %[~,out] = system(['fdp -Tpdf ' out_file ' -o ' out_fig '_fdp.pdf'])
    %[~,out] = system(['neato -Tpdf ' out_file ' -o ' out_fig '_neato.pdf'])
    [~,out] = system(['dot -Tpdf ' out_file ' -o ' out_fig '_dot.pdf'])
    %[~,out] = system(['twopi -Tpdf ' out_file ' -o ' out_fig '_twopi.pdf'])
    [~,out] = system(['circo -Tpdf ' out_file ' -o ' out_fig '_circo.pdf'])
end

