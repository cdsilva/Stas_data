function F_l_new = rotate_fhat(F_l, lmax, angles)
                
alpha = angles(1);
beta = angles(2);
gamma = angles(3);

D = zeros(length(F_l));
max_ind = 1;
for l=0:lmax
    D(max_ind:max_ind+2*l, max_ind:max_ind+2*l) = wigner_big_D(l, alpha, beta, gamma);
    max_ind = max_ind+2*l+1;
end
F_l_new = F_l * D;