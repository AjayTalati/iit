function perspective = comp_pers_cpt(num_nodes_indices,denom_nodes_indices,numerator_state,bf_option)

%  compute BRs and FRs for a single perspective but given some fixed
%  current state

if isempty(denom_nodes_indices)
    perspective = [];
    return
% elseif isempty(num_nodes_indices)
% %     num_sys_nodes = denom_nodes_indices(1).num_sys_nodes;
% %     denom_conditional_joint_size = ones(1,2*num_sys_nodes);
% %     denom_conditional_joint_size(1:num_sys_nodes == denom_nodes_indices
%     denom_conditional_joint = [];
%     return
end

global nodes

num_sys_nodes = nodes(1).num_sys_nodes;

if strcmp(bf_option,'backward')
    
    denom_nodes = nodes(denom_nodes_indices);
    num_nodes_shift = num_nodes_indices + num_sys_nodes;
    numerator_nodes = nodes(num_nodes_shift);
    
    % no nodes in numerator means maxent over denom
    if isempty(num_nodes_indices)
        
        perspective_dim_sizes = ones(1,num_sys_nodes);
        perspective_dim_sizes(denom_nodes_indices) = [denom_nodes.num_states];
        perspective = ones(perspective_dim_sizes)./prod(perspective_dim_sizes);
        return
        
    end
    
    % setup cell array for conditioning
    conditioning_indices = cell(1,2*num_sys_nodes);
    for i = 1:2*num_sys_nodes
        conditioning_indices{i} = ':';
    end
    
    % we do the first iteration outside the for main foor loop so we can
    % initialize the joint
    this_node_conditioning_indices = conditioning_indices;
    this_node_conditioning_indices{numerator_nodes(1).num} = numerator_state(numerator_nodes(1).num - num_sys_nodes) + 1;
    numerator_conditional_joint = numerator_nodes(1).cpt(this_node_conditioning_indices{:});
    
    prob_current_state = sum(numerator_conditional_joint(:));
    
    % marginalize over nodes not in denominator, these nodes are outside the
    % system for this iteration or they are outside a partition - either
    % way we apply maxent prior + marginalization
    for j = 1:num_sys_nodes
        
        if ~any(j == denom_nodes_indices) && any(j == numerator_nodes(1).input_nodes)
            numerator_conditional_joint = ...
                sum(numerator_conditional_joint,j)./size(numerator_conditional_joint,j);
        end
        
    end
    
    for i = 2:length(num_nodes_indices)
        
        this_node_conditioning_indices = conditioning_indices;
        this_node_conditioning_indices{numerator_nodes(i).num} = numerator_state(numerator_nodes(i).num - num_sys_nodes) + 1;
        next_num_node_distribution = numerator_nodes(i).cpt(this_node_conditioning_indices{:});
        
        prob_current_state = prob_current_state * sum(next_num_node_distribution(:));

        % marginalize over nodes not in denom, these nodes are outside the
        % system for this iteration or they are outside a partition - either
        % way we apply maxent prior/marginalization
        for j = 1:num_sys_nodes

            if ~any(j == denom_nodes_indices) && any(j == numerator_nodes(i).input_nodes)
                next_num_node_distribution = ...
                    sum(next_num_node_distribution,j)./size(next_num_node_distribution,j);
            end
        end
        
        % the magic
        numerator_conditional_joint = bsxfun(@times,numerator_conditional_joint,next_num_node_distribution);
    end
    
    % conditioning on fixed nodes
%     perspective = numerator_conditional_joint .* (1/prod([denom_nodes.num_states])) ./ prob_current_state;
    perspective = numerator_conditional_joint ./ sum(numerator_conditional_joint(:));
    
    
    
% P(denom_nodes_f | num_nodes_c = numerator_state) = P(denom_nodes_c | num_nodes_p = numerator_state)
elseif strcmp(bf_option,'forward')
    
    denom_nodes_shift = denom_nodes_indices + num_sys_nodes;
    denom_nodes = nodes(denom_nodes_shift);
    
% % % %     % no nodes in numerator means maxent over denom
% % % %     if isempty(num_nodes_indices)
% % % %         
% % % %         perspective_dim_sizes = ones(1,num_sys_nodes);
% % % %         perspective_dim_sizes(denom_nodes_indices) = [denom_nodes.num_states];
% % % %         perspective = ones(perspective_dim_sizes)./prod(perspective_dim_sizes);
% % % %         return
% % % %         
% % % %     end
% % % %     
    % we do the first iteration outside the for main foor loop so we can
    % initialize the joint
    denom_conditional_joint = denom_nodes(1).cpt;
    
    conditioning_indices = cell(1,2*num_sys_nodes);
    for i = 1:2*num_sys_nodes
        conditioning_indices{i} = ':';
    end
    
    % marginalize over nodes not in numerator, these nodes are outside the
    % system for this iteration or they are outside a partition - either
    % way we apply maxent prior/marginalization
    for j = 1:num_sys_nodes
        
        if ~any(j == num_nodes_indices) && any(j == denom_nodes(1).input_nodes)
            denom_conditional_joint = ...
                sum(denom_conditional_joint,j)./size(denom_conditional_joint,j);
        elseif any(j == num_nodes_indices)
            conditioning_indices{j} = numerator_state(j) + 1;
        end
        
    end
    
    for i = 2:length(denom_nodes)
        
        next_denom_node_distribution = denom_nodes(i).cpt;
        
        % marginalize over nodes not in denom, these nodes are outside the
        % system for this iteration or they are outside a partition - either
        % way we apply maxent prior/marginalization
        for j = 1:num_sys_nodes

            if ~any(j == num_nodes_indices) && any(j == denom_nodes(i).input_nodes)
                next_denom_node_distribution = ...
                    sum(next_denom_node_distribution,j)./size(next_denom_node_distribution,j);
            end
        end
        
        % the magic
        denom_conditional_joint = bsxfun(@times,denom_conditional_joint,next_denom_node_distribution);
    end
    
    
    % conditioning on fixed nodes
    denom_conditional_joint = denom_conditional_joint(conditioning_indices{:});
    permute_order = [num_sys_nodes+1:2*num_sys_nodes 1:num_sys_nodes];
    perspective = permute(denom_conditional_joint,permute_order);

end