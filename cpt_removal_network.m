function network = cpt_removal_network(network)
    % subset - build a cell array that contains all of the subsets
    % subsets builds arrays that use the actual node numbers as opposed to
    % logicals - perhaps we should make one of these that is global as well
    J = network.connect_mat;
    subsys_conn_mat = J(this_subset, this_subset);


    %Larissa: New tpm, new logic types! and new nodes
    %Maybe all I need is the cpt! -> No rather everything to be save

    % setup node strucs and cpts for each node
    %inputs = struct('num',{1 2},'name',{'A_p' 'B_p'},'num_states',{2 2},'state_names',{{'0' '1'}},'logic_type',{2 3})
    logic_types = get(handles.logic_types,'Data');
    % init struct array
    nodes(2*num_nodes) = struct('num',2*num_nodes,'name',[num2str(num_nodes) '_c'],'num_states',2,...
                                'state_names',{{'0' '1'}},'logic_type',logic_types(num_nodes),'cpt',[],...
                                'num_sys_nodes',num_nodes,'input_nodes',[]);
    % make past node structs                        
    for i = 1:num_nodes
        nodes(i) = struct('num',i,'name',[num2str(i) '_p'],'num_states',2,...
                                'state_names',{{'0' '1'}},'logic_type',logic_types(i),'cpt',[],...
                                'num_sys_nodes',num_nodes,'input_nodes',[]);
    end

    % make current node structs and their tpms
    for i = 1:num_nodes
        nodes(num_nodes + i) = struct('num',num_nodes + i,'name',[num2str(i) '_c'],'num_states',2,...
                                'state_names',{{'0' '1'}},'logic_type',logic_types(i),'cpt',[],...
                                'num_sys_nodes',num_nodes,'input_nodes',[]);
        input_nodes = 1:num_nodes;
        input_nodes_indices = input_nodes(logical(connectivity_matrix(i,:)));
        nodes(num_nodes + i).input_nodes = input_nodes_indices;

    %     input_nodes = nodes(input_nodes_indices);
    %     nodes(num_nodes + i).cpt = cpt_factory_mechs(nodes(num_nodes + i),input_nodes,2*num_nodes,noise);
    %     disp(nodes(num_nodes + i).cpt)
    %     test_cpt = cpt_factory_tpm(nodes(num_nodes + i), input_nodes_indices, nodes, 2*num_nodes, tpm);
        nodes(num_nodes + i).cpt = cpt_factory_tpm(nodes(num_nodes + i), input_nodes_indices, nodes, 2*num_nodes, tpm);

    %     if any(nodes(num_nodes + i).cpt ~= test_cpt)
    %         disp('error')
    %     end



    end
end
