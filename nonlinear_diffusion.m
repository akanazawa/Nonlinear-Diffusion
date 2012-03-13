%%%%%%%%%%
% Non-linear diffusion method in "A review of
% nonlinear diffusion filtering," by Joachim Weickert. eq (27),
% (28), with eq(3)
% Construct D with eigenvectors lam1 = g(|\grad u|^2), lam2 = 1
% do D_tu =  div(D.* \grad u)
%%%%%%%%%%
function [u] = nonlinear_diffusion(I, lam, sig)
    itr = 200;
    diffusivity = .1;
    u = I;
    lam2 = 1;
    %    filter = fspecial('gaussian', ceil(sig)*6, sig);
    filter = fspecial('gaussian', ceil(sig)*3, sig);
    for t = 1:itr
        [Ix Iy] = gradient(u);
        [Ix_sig Iy_sig] = gradient(imfilter(u,filter, ...
                                            'replicate'));
        Du_x = zeros(size(I));
        Du_y = zeros(size(I));
        normGrad = Ix_sig.^2 + Iy_sig.^2;
        g = 1./(1 +(normGrad/lam^2));
        for s = 1:numel(I) % do it pixel wise
            gradU = [Ix(s); Iy(s)];
            gradU_sig = [Ix_sig(s) ; Iy_sig(s)];            
            % construct D. v1 is parallell to grad u, v2 is orth
            v1 = gradU_sig./norm(gradU_sig);
            v2 = [v1(2) ; -v1(1)];
            D = [v1 v2]*[g(s) 0; 0 1]*[v1 v2]';
            % D = g(s)*v1*v1' + 1*v2*v2';
            Du = D*gradU;  
            Du_x(s) = Du(1);
            Du_y(s) = Du(2);            
        end               
        ut = divergence(Du_x, Du_y);
        u = u + diffusivity*ut;        
    end
end
