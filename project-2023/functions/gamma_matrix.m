%% function to generate a gamma matrix
function [matrix] = gamma_matrix(p_x, K, k_hat, tau, p_z)

% define k indices from -(K-1)/2 to (K-1)/2
k_s = linspace(-(K-1)/2, (K-1)/2, K);
ck_s = (4/3) * tau * k_s;

% constant term in front of the vector
alpha = k_hat * exp(-p_z * (pi/tau));

% create gamma matrix
matrix = zeros(3, K);
for i = 1:K
    c_k = ck_s(i);
    const = (pi/tau) * (p_x - c_k);
    matrix(1, i) = sin(const);
    matrix(2, i) = cos(const);
    matrix(3, i) = p_z * sin(const) - (p_x - c_k) * cos(const);
end
% multiply by alpha
matrix = alpha * matrix;
end