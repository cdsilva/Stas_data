function [f_new, angles_opt] = align_fhat(f, f_template, lmax)

alpha_range = -pi/8:pi/16:pi/8;
beta_range = 0:pi/8:(2*pi-pi/8);
gamma_range = -pi/8:pi/16:pi/8;

angles0 = zeros(1,3);
err = inf;
for a=alpha_range
    for b=beta_range
        for c=gamma_range
            tmp_err = sse_f(f_template, f, lmax, [a b c]);
            if tmp_err < err
                err = tmp_err;
                angles0 = [a b c];
            end
        end
    end
end
        
angles_opt = fminsearch(@(angles) sse_f(f_template, f, lmax, angles), angles0);

f_new = rotate_fhat(f, lmax, angles_opt);


function y = sse_f(f_template, f, lmax, angles)

y = norm(f_template - rotate_fhat(f, lmax, angles));