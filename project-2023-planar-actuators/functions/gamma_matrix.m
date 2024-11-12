%% function to generate a gamma matrix
function [matrix] = gamma_matrix(p_x, K, k_hat, tau, p_z)

% define k indices from -(K-1)/2 to (K-1)/2
k_s = linspace(-(K-1)/2, (K-1)/2, K);
ck_s = (4/3) * tau * k_s;

% constant term in front of the vector
const = k_hat * exp(-p_z * (pi/tau));

% create gamma matrix
matrix = zeros(3, K);
for i = 1:K
    c_k = ck_s(i);
    alpha_k = (pi/tau) * (p_x - c_k);
    matrix(1, i) = sin(alpha_k);
    matrix(2, i) = cos(alpha_k);
    matrix(3, i) = p_z * sin(alpha_k) - (p_x - c_k) * cos(alpha_k);
end
% multiply by const
matrix = const * matrix;
end