% Solve the minmax abs current problem
%
% inputs:
%   - Gamma: 3xK Gamma coupling matrix
%   - w_des: 3x1 wrench vector with desired forces
function i_com = min_max_abs_cur(Gamma, w_des)
    options = optimoptions('linprog','Display','none');

    % Get number of coils
    K = length(Gamma);

    % add one row of zeros to the eq constraint for variable x(K+1) = M
    % so we have: [i_1, i_2, ..., i_k, M]
    Gm_px = [Gamma zeros([length(w_des),1])];
    
    % equality constraint
    Aeq = Gm_px;
    beq = w_des;
    
    % ineq constraint 
    % variable x(K+1) = M
    % -1*i_k - M <= 0
    neg_ones = -1*ones([1,K+1]);
    Aineq_neg = diag(neg_ones);
    Aineq_neg(:, end) = -1;
    Aineq_neg = Aineq_neg(1:end-1,:);
    
    % i_k - M <= 0
    pos_ones = -1*neg_ones;
    Aineq_pos = diag(pos_ones);
    Aineq_pos(:, end) = -1;
    Aineq_pos = Aineq_pos(1:end-1,:);
    
    % merge Aineq_pos and Aineq_neg
    Aineq = [Aineq_pos; 
             Aineq_neg];
    bineq = zeros(2*K,1);
    
    % define the objective = M = x(K+1)
    f = zeros([1,K+1]);
    f(K+1) = 1;
    
    % solve for i
    [x,fval,exitflag,output,lambda] = linprog(f, Aineq, bineq, Aeq, beq, [], [], options);
    i_com = x(1:K);
end