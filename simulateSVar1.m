function [SV, h] = simulateSVar1(hrho, hbar, hT, hvol, Nfedraws, Nhorizons,rndStream)
N        = length(hrho);
h        = NaN(N, Nfedraws, Nhorizons);
hshocks  = hvol * randn(rndStream, N, Nfedraws .* Nhorizons);
hshocks  = reshape(hshocks, N, Nfedraws, Nhorizons);
intercept = (1 - hrho) .* hbar;
j = 1;
h(:,:,j) = intercept + hrho .* hT + hshocks(:,:,j);
for j = 2 : Nhorizons
    h(:,:,j) = intercept + hrho .* h(:,:,j-1) + hshocks(:,:,j);
end
SV = exp(h * .5);