function Make_Graph_for_graphviz(Node,Edge,Graph,out_file,out_fig,my_mot)

    disp('Plotting Graph')

    %brain_layer_marker = textread('brain_layer_marker.txt','%s','delimiter','\n');

    font_size = '50';
    linewidth = '14';
    arrow_size = '3';
    min_len = '1';

    N_Node = size(Node,1);
    N_Edge = size(Edge,1);

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
		label = num2str(Edge.ll(i),2);
        %if Node.TF{tar_ind}.Count > 1
            %label = Edge.gene{i};
        %else
            %label = '';
        %end
        fprintf(fid,'\t%s -> %s [label="%s",arrowsize="%s",arrowhead="%s",fontsize=%s,minlen=%s];\n',Edge.source{i},Edge.target{i},label,arrow_size,arrow_head,font_size,min_len);
    end

    % Create Nodes
    for i = 1:N_Node

        % Make html table for brain layer markers
        %label = ['<<table border="0" cellspacing="0"><tr><td>' Node.name{i} '</td></tr>'...
        %         '<tr><td><IMG SRC="' Node.image{i} '" STYLE="width:75px;height:50px;"/></td></tr>'];
        %label = [label '<tr><td>' Node.TF{i} '</td></tr>'];
        %label = [label '</table>>'];
        if strcmp(my_mot,'Genes')
            Node.image{i} = '';
        end
        if strcmp(my_mot,'Genes')
            label = Node.name{i};     
        else
            label = ['<<table border="0" cellspacing="0"><tr><td>' Node.name{i} '</td></tr>'];
            if Node.TF{i}.Count > 0
                label = [label '<tr><td>'];
                for tf = Node.TF{i}.keys
                    label = [label tf{:} ':' num2str(Node.TF{i}(tf{:})) ' | '];
                end
                label(end-2:end) = '';
                label = [label '</td></tr>'];
            end
            for go = Node.go{i}
                label = [label '<tr><td>' go{:} '</td></tr>'];
            end
            label = [label '</table>>'];
        end
        
		
		
		
		%fprintf(fid,'\t%s [shape=%s,label=%s,image="%s",color=%s,style="setlinewidth(%s)",fontsize=20,style="filled", fillcolor=%s];\n',Node.name{i},Node.shape{i},label,Node.image{i},Node.color{i},linewidth,Node.color{i}); 
        fprintf(fid,'\t%s [shape=%s,label=%s,image="%s",color=%s,style="setlinewidth(%s)",fontsize=%s];\n',Node.name{i},Node.shape{i},label,Node.image{i},Node.color{i},linewidth,font_size);
        %fprintf(fid,'\t%s [shape=%s,label=%s,image="%s",color=%s,style="setlinewidth(%s)",style="filled",fillcolor="0 0 .5"];\n',Node.name{i},Node.shape{i},label,Node.image{i},Node.color{i},linewidth);
    end
    fprintf(fid,'}\n');

    fclose(fid);

    %[~,out] = system(['fdp -Tpdf ' out_file ' -o ' out_fig '_fdp.pdf'])
    %[~,out] = system(['neato -Tpdf ' out_file ' -o ' out_fig '_neato.pdf'])
    [~,out] = system(['dot -Tpdf ' out_file ' -o ' out_fig '_dot.pdf'])
    %[~,out] = system(['twopi -Tpdf ' out_file ' -o ' out_fig '_twopi.pdf'])
    [~,out] = system(['circo -Tpdf ' out_file ' -o ' out_fig '_circo.pdf'])
    [~,out] = system(['circo -Tsvg ' out_file ' -o ' out_fig '_circo.svg'])
end

