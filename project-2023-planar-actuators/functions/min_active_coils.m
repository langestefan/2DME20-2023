% Solve the min number of active coils problem
%
% inputs:
%   - Gamma: 3xK Gamma coupling matrix
%   - w_des: 3x1 wrench vector with desired forces

function i_com = min_active_coils(Gamma, w_des)
    options = optimoptions('linprog','Display','none');

    % Get number of coils
    K = length(Gamma);

    % inequality constraints
    bineq = zeros(2*K,1);
    diag_nones = diag(-1 * ones(1,K));
    Aineq = [diag_nones   zeros([K,K]);
             zeros([K,K]) diag_nones];
    
    % equality constraints
    beq = w_des;
    
    % objective f = sum(i_plus) + sum(i_neg)
    f = ones([2*K,1]);

    % update equality constraints
    Aeq = [Gamma -Gamma];

    % solve for i
    [x,fval,exitflag,output,lambda] = linprog(f, Aineq, bineq, Aeq, beq, [], [], options);

    % transform i+ and i- back to the original i
    % i = i_pos - i_neg
    i_com = x(1:K) - x(K+1:end);    
end