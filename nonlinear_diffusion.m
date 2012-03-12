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
    filter = fspecial('gaussian', ceil(sig)*6, sig);
    for t = 1:itr
        [Ix Iy] = gradient(u);
        [Ix_sig Iy_sig] = gradient(imfilter(u,filter,'replicate'));
        Du_x = zeros(size(I));
        Du_y = zeros(size(I));
        for s = 1:numel(I) % do it pixel wise
            gradU = [Ix(s) ; Iy(s)];
            gradU_sig = [Ix_sig(s) ; Iy_sig(s)];            
            % construct D. v1 is parallell to grad u, v2 is orth
            %http://campar.in.tum.de/twiki/pub/Chair/TeachingSs04ImageSeg/VL_Segmentierung_kap8.pdf
            if norm(gradU_sig) ~= 0
                lam1 = 1./(1+(norm(gradU_sig)/lam).^2);
                v1 = gradU_sig/norm(gradU_sig);
                v2 = [v1(2) ; -v1(1)];
                D = [v1 v2]*[lam1 0; 0 1]*[v1 v2]';
                %D = lam1*v1*v1' + lam2*v2*v2';
                Du = D*gradU;                
                Du_x(s) = Du(1);
                Du_y(s) = Du(2);
            end
        end
        ut = divergence(Du_x, Du_y);
        u = u + diffusivity*ut;        
    end
end
