%%%%%%%%%%
% Regularized Perona-Malik diffusion method from Weickert paper "A review of
% nonlinear diffusion filtering," by Joachim Weickert. eq (23),
% (21)
% Du/Dt = div(g(|grad(u_sig)|^2 .* grad u)
%%%%%%%%%%
function [u] = perona_malik(I, lam, sig)
    itr = 200;
    diffusivity = .1;
    u = I;
    h = fspecial('gaussian', ceil(sig)*5, sig);
    for t = 1:itr
        [Ix Iy] = gradient(u);
        [Ix_sig Iy_sig] = gradient(imfilter(u,h, 'replicate'));
        normGrad = Ix_sig.^2 + Iy_sig.^2;
        g = 1./(1+ (normGrad/lam^2));
        ut = divergence(g.*Ix, g.*Iy);
        u = u + diffusivity*ut;        
    end
end
