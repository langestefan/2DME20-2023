% Solve the minimum absolute power minimization problem
%
% inputs:
%   - R: KxK diagonal resistance matrix
%   - Gamma: 3xK Gamma coupling matrix
%   - w_des: 3x1 wrench vector with desired forces
function i_com = min_abs_power(R, Gamma, w_des)
    % Get number of coils
    K = length(Gamma);

    % need extra zeros for the equation
    z_vec = zeros(K, 1);
    z_mtx = zeros(K-16, K-16);
    x = [z_vec; w_des];

    % obtain optimal commutation current i_com for w_desired
    A = [R      Gamma' ;
         Gamma  z_mtx ];
    
    % solve Ax = b 
    % where x = [i; mu] with i = (Kx1) and mu = (3x1)
    i_com = A\x;
    i_com = i_com(1:K); % take the i_com part of x
end