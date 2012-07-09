function [phi_MIP prob prob_prod_MIP MIP] = phi_comp_bORf(options,x0,x,p,bf,b_table,x0_s)
% Larissa: for smart purviews, op_context is assumed 0, op_min is assumed
% bf is back/forward flag (back = 1, forward = 2)
op_fb = options(1);
op_phi = options(2);
op_disp = options(3);
op_single = options(4);
op_ex = options(5);
op_context = options(6);
op_whole = options(7);
op_min = options(9);
op_normalize = options(14);
op_small_phi = options(16);

global BRs, global FRs


N = length(x);
N0 = length(x0);
%% unpartitioned transition repertoire
% prob_w_old = Rs{convi(x0),convi(x)};

    current = convi(x0); other = convi(x);
    
    if (bf == 1)
        if isempty(BRs{current,other})
            BRs{current,other} = comp_pers_single(current,other,x0_s,p,b_table,bf);
        end
        prob_w = BRs{current,other};
    elseif (bf == 2)
        if isempty(FRs{current,other})
            FRs{current,other} = comp_pers_single(current,other,x0_s,p,b_table,bf);
        end
        prob_w = FRs{current,other};
    end
    
%     disp('CHECK NEW COMPUTATION:');
%     disp('prob_w_old:')
%     disp(prob_w_old);
%     disp('prob_w:')
%     disp(prob_w);
%     disp(prob_w_old == prob_w);

    
prob = cell(2,1);
prob{bf} = prob_w;

%% more than one
if N ~= 0
    [x_b1 x_b2 N_b] = bipartition(x,N); % partition of xp
else
    x_b1{1} = []; x_b2{1} = []; N_b = 1;
end

[x0_b1 x0_b2 N0_b] = bipartition(x0,N0,1); % partition of x0

phi_cand = zeros(N_b,N0_b,2,2);
prob_prod_vec = cell(N_b,N0_b,2,2);

for i=1: N_b % past or future
    x_1 = x_b1{i};
    x_2 = x_b2{i};
    for j=1: N0_b % present
        x0_1 = x0_b1{j};
        x0_2 = x0_b2{j};
        Norm = Normalization(x_1,x_2,x0_1,x0_2);

        if Norm ~= 0
            current_1 = convi(x0_1); current_2 = convi(x0_2);
            other_1 = convi(x_1); other_2 = convi(x_2);
            
            if (bf == 1)
                if isempty(BRs{current_1,other_1})
                    BRs{current_1,other_1} = comp_pers_single(current_1,other_1,x0_s,p,b_table,bf);
                end
                prob_p1 = BRs{current_1,other_1};
                
                if isempty(BRs{current_2,other_2})
                    BRs{current_2,other_2} = comp_pers_single(current_2,other_2,x0_s,p,b_table,bf);
                end
                prob_p2 = BRs{current_2,other_2};
                
            elseif (bf == 2)
                if isempty(FRs{current_1,other_1})
                    FRs{current_1,other_1} = comp_pers_single(current_1,other_1,x0_s,p,b_table,bf);
                end
                prob_p1 = FRs{current_1,other_1};
                
                if isempty(FRs{current_2,other_2})
                    FRs{current_2,other_2} = comp_pers_single(current_2,other_2,x0_s,p,b_table,bf);
                end
                prob_p2 = FRs{current_2,other_2};
            end
            
            prob_p = prob_prod_comp(prob_p1,prob_p2,x,x_1,0);
            if (op_small_phi == 0)
                phi = KLD(prob{bf},prob_p);
            elseif (op_small_phi == 1)
                phi = emd_hat_gd_metric_mex(prob{bf}',prob_p',gen_dist_matrix(length(prob_p)));
            end
            prob_prod_vec{i,j,bf} = prob_p;
            
        else
            prob_prod_vec{i,j,bf} = [];
            phi = Inf;
        end

        phi_cand(i,j,bf,1) = phi;
        phi_cand(i,j,bf,2) = phi/Norm;

            % fprintf('phi=%f phi_norm=%f %s-%s -%s\n',phi,phi/Norm,mod_mat2str(xp_1),mod_mat2str(x0_1),mod_mat2str(xf_1));
    end
end

MIP = cell(2,2,2);
phi_MIP = zeros(1,2);
prob_prod_MIP = cell(2,1);

[phi_MIP(bf) i j] = min2(phi_cand(:,:,bf,1),phi_cand(:,:,bf,2),op_normalize);
prob_prod_MIP{bf} = prob_prod_vec{i,j,bf};

MIP{1,1,bf} = x_b1{i};
MIP{2,1,bf} = x_b2{i};
MIP{1,2,bf} = x0_b1{j};
MIP{2,2,bf} = x0_b2{j};

end

function Norm = Normalization(xp_1,xp_2,x0_1,x0_2,xf_1,xf_2)

if nargin == 4
    Norm = min(length(x0_1),length(xp_2)) + min(length(x0_2),length(xp_1));
else
    Norm = min(length(x0_1),length(xp_2)) + min(length(x0_2),length(xp_1)) ...
        + min(length(x0_1),length(xf_2)) + min(length(x0_2),length(xf_1));
end

end

function [X_min i_min j_min k_min] = min3(X,X2)
X_min = Inf; % minimum of normalized phi
X_min2 = Inf; % minimum of phi
i_min = 1;
j_min = 1;
k_min = 1;

for i=1: size(X,1)
    for j=1: size(X,2)
        for k=1: size(X,3)
            if X(i,j,k) <= X_min && X2(i,j,k) <= X_min2
                X_min = X(i,j,k);
                X_min2 = X2(i,j,k);
                i_min = i;
                j_min = j;
                k_min = k;
            end            
        end
    end
end

end


function [phi_min_choice i_min j_min] = min2(phi,phi_norm,op_normalize)
phi_norm_min = Inf; % minimum of normalized phi
phi_min = Inf; % minimum of phi
i_min = 1;
j_min = 1;

if (op_normalize == 1 || op_normalize == 2)
    for i=1: size(phi,1)
        for j=1: size(phi,2)
%             if phi_norm(i,j) <= phi_norm_min && phi(i,j) <= phi_min
            if phi_norm(i,j) <= phi_norm_min
                phi_min = phi(i,j);
                phi_norm_min = phi_norm(i,j);
                i_min = i;
                j_min = j;
            end
        end
    end
else
    for i=1: size(phi,1)
        for j=1: size(phi,2)
            if phi(i,j) <= phi_min
                phi_min = phi(i,j);
                phi_norm_min = phi_norm(i,j);
                i_min = i;
                j_min = j;
            end
        end
    end
end

if (op_normalize == 0 || op_normalize == 1)
    phi_min_choice = phi_min;
else
    phi_min_choice = phi_norm_min;
end

end
