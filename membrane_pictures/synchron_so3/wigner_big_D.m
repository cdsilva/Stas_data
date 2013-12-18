function D = wigner_big_D(l, alpha, beta, gamma)

D = zeros(2*l+1);

fm = 0:2*l;
fm(1) = 1;
fm = cumprod(fm);

for m=-l:l
    for n=-l:l
        s = max([0, m-n]):min([l+m, l-n]);
        
        d = sqrt(fm(l+m+1) * fm(l-m+1) * fm(l+n+1) * fm(l-n+1)) * sum((-1).^(n-m+s)./(fm(l+m-s+1).*fm(s+1).*fm(n-m+s+1).*fm(l-n-s+1)).*cos(beta/2).^(2*l+m-n-2*s).*sin(beta/2).^(n-m+2*s));
        D(m+l+1,n+l+1) = exp(-sqrt(-1)*m*alpha) * d *exp(-sqrt(-1)*n*gamma);
    end
end


