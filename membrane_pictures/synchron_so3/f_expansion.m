function f1 = f_expansion(theta, phi, F_l, l_max)

F_l = reshape(F_l, 1, []);

f1 = zeros(size(theta));
for l=0:l_max
    start_idx = sum(2*(0:(l-1))+1)+1;
    sph_harm = Y_lm(l, reshape(theta, 1, []), reshape(phi, 1, []));
    f1 = f1 + reshape(sum(repmat(F_l(start_idx:start_idx+2*l).', 1, size(sph_harm,2)).*conj(sph_harm), 1), size(f1));
end