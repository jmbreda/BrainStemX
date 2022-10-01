function Make_Graph_for_graphviz(Node,Edge,Graph,out_file,out_fig,my_mot)

    disp('Plotting Graph')

    %brain_layer_marker = textread('brain_layer_marker.txt','%s','delimiter','\n');

    font_size = '32';
    linewidth = '8';
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
		label = num2str(Edge.ll(i),2);
        arrow_head = 'normal';
        tar_ind = find(strcmp(Node.name,Edge.target{i}));
        if Node.TF{tar_ind}.Count > 0
            node_corr = Node.TF{tar_ind}(Edge.gene{i});
            if sign(node_corr)==-1
                arrow_head = 'inv';
            end
        end
        %fprintf(fid,'\t%s -> %s [arrowsize="%s",arrowhead="%s",fontsize=%s,minlen=%s,headlabel="%s",taillabel="%s"];\n',Edge.source{i},Edge.target{i},arrow_size,arrow_head,font_size,min_len,Edge.source{i},Edge.gene{i});
        fprintf(fid,'\t%s -> %s [label="%s",arrowsize="%s",arrowhead="%s"];\n',Edge.source{i},Edge.target{i},label,arrow_size,arrow_head);
    end


    % Create Nodes
    if strcmp(my_mot,'Genes') || strcmp(my_mot,'All')
        All_Nodes = Node;

        fprintf(fid,'\tsubgraph cluster%u {\n',0);
        fprintf(fid,'\t\tstyle=invis;\n');
        ind_node = find(All_Nodes.phase==2);
        ind_node = intersect(ind_node,find(All_Nodes.in_degree==0));
        Node = All_Nodes(ind_node,:);

        for i = 1:size(Node,1);

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
        fprintf(fid,'\t}\n');


        Phases = {'Expansion','Neurogenesis','Gliogenesis'}; 

        for p = 1:3
            % New subgragh
            fprintf(fid,'\tsubgraph cluster%u {\n',p);
            fprintf(fid,'\t\tlabel = "%s" fontsize=%s style=invis;\n',Phases{p},font_size);

            ind_node = find(All_Nodes.phase==p);
            if p==2
                ind_node = intersect(ind_node,find(All_Nodes.in_degree>0));
            end

            Node = All_Nodes(ind_node,:);

            for i = 1:size(Node,1);

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
            fprintf(fid,'\t}\n');
        end
    else
        for i = 1:size(Node,1);

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
                for go = Node.go{i}
                    label = [label '<tr><td>' go{:} '</td></tr>'];
                end
                label = [label '</table>>'];
            end 
            fprintf(fid,'\t\t%s [shape=%s,label=%s,image="%s",color=%s,style="setlinewidth(%s)",fontsize=%s];  \n',Node.name{i},Node.shape{i},label,Node.image{i},Node.color{i},linewidth,font_size);
        end
    end
    fprintf(fid,'}\n');

    fclose(fid);

    %[~,out] = system(['fdp -Tpdf ' out_file ' -o ' out_fig '_fdp.pdf'])
    %[~,out] = system(['neato -Tpdf ' out_file ' -o ' out_fig '_neato.pdf'])
    [~,out] = system(['dot -Tpdf ' out_file ' -o ' out_fig '_dot.pdf'])
    [~,out] = system(['dot -Tsvg ' out_file ' -o ' out_fig '_dot.svg'])
    %[~,out] = system(['twopi -Tpdf ' out_file ' -o ' out_fig '_twopi.pdf'])
    [~,out] = system(['circo -Tpdf ' out_file ' -o ' out_fig '_circo.pdf'])
    [~,out] = system(['circo -Tsvg ' out_file ' -o ' out_fig '_circo.svg'])
end

