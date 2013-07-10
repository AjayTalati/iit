function TPMshaped = cpt_factory_tpm_SxS(tpm, El_num_states, num_sys_nodes)

TPMshaped.FR = reshape(tpm, [El_num_states, prod(El_num_states)]);
%TPMshaped.BR = reshape(tpm);


% 
% dim_sizes = ones(1,num_sys_nodes);
% dim_sizes(1:num_sys_nodes) = El_num_states;
% cptTPM = cell(prod(El_num_states),1);
% 
% for s = 1:size(tpm,1)
%     prob_dist = reshape(tpm(s,:),dim_sizes);
%     cptTPM{s} = prob_dist;
% end
% cptTPM = reshape(cptTPM, El_num_states);
end