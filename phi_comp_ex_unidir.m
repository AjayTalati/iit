function [max_phi_MIP, prob, j_max, network] = phi_comp_ex_unidir(subsystem,M1,M2,numerator,whole_sys_state,network,bf_option,bfcut_option)
%pf_tag is 1 for past and 2 for future

op_extNodes = network.options(11);

nodes_vec = subsystem;
N = numel(nodes_vec);
num_states_subsys = prod([network.nodes([subsystem]).num_states]);
subsets_subsys = cell(num_states_subsys-1,1); % subtract one since we don't consider the empty system
for i = 1:num_states_subsys-1 % don't include empty set, this is why for-loop starts at 2
    subsets_subsys{i} = nodes_vec(logical(network.b_table{i+1,N}));
end

if op_extNodes == 0
    extNodes = setdiff(network.full_system, subsystem);
else
    extNodes = [];
end 

phi_MIP = zeros(num_states_subsys-1,2);
prob_cand = cell(num_states_subsys-1,1);

for i=1: num_states_subsys-1
    %Larissa smart purviews: Only test those connections that actually exist
    % could made even smarter here for unidirectional noise
    denom = subsets_subsys{i};
    if strcmp(bf_option,'backward')
        if nnz(sum(network.connect_mat(numerator,denom),1) == 0) == 0 % all denom is input of numerator (numerator) --> phiBR
            [phi_temp, prob_temp ,~,~, network] = phi_comp_bORf(subsystem,numerator,denom,whole_sys_state,network,1,M1,M2,bfcut_option);
            phi_MIP(i) = phi_temp(1);
            prob_cand{i} = prob_temp{1};
        else 
            uniform_dist = ones(num_states_subsys,1)/num_states_subsys;
            prob_cand{i} = uniform_dist;
        end
    elseif strcmp(bf_option,'forward')
        if nnz(sum(network.connect_mat(denom,numerator),2) == 0) == 0 % denom is output
            [phi_temp, prob_temp ,~,~, network] = phi_comp_bORf(subsystem,numerator,denom,whole_sys_state,network,2,M1,M2,bfcut_option);
            phi_MIP(i) = phi_temp(2);
            prob_cand{i} = prob_temp{2};
        else 
            forward_maxent_dist = comp_pers_cpt(network.nodes,[],subsystem,whole_sys_state,'forward',extNodes);
            prob_cand{i} = forward_maxent_dist(:);
        end
    end  
end

% now take the largest phi, if equal take the bigger set
[max_phi_MIP j_max] = max_ex(phi_MIP,subsets_subsys);
denom = subsets_subsys{j_max};
prob = prob_cand{j_max};
if length(denom) ~= N
    if strcmp(bf_option,'backward')
        prob = expand_prob(prob,subsystem,denom);
    elseif strcmp(bf_option,'forward')
        denom_rest = pick_rest(subsystem,denom);
        fmaxent_denom_rest = comp_pers_cpt(network.nodes,[],denom_rest,whole_sys_state,'forward',extNodes);
        prob = expand_prob_general(prob,subsystem,denom,fmaxent_denom_rest(:));
    end
end
end

%% subfunction
function [X_max i_max] = max_ex(X,subsets_subsys)
% exclusion principle: if the value is the same, take the bigger one
epsilon = 10^-6;    %EMD has very low precision
X_max = -Inf;
i_max = 1;
s_max = 0;
for i=1: size(X,1)
    s = length(subsets_subsys{i});
    cond1 = X(i) > X_max;
    cond2 = abs(X(i) - X_max) < epsilon && s>= s_max;
    if cond1 || cond2
        X_max = X(i);
        i_max = i;
        s_max = s;
    end
end

end