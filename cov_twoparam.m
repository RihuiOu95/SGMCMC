function [res] = cov_twoparam(phi1, phi2, rho)
res = mvnpdf([phi1,phi2],[0,0],[1 rho;rho 1]) * normpdf(phi1) .* normcdf((phi2 - rho.*phi1)./sqrt(1-rho^2)) ./ mvncdf([phi1,phi2],[0,0],[1 rho;rho 1]) ...
    - mvnpdf([phi1,-phi2],[0,0],[1 -rho;-rho 1]) * normpdf(phi1) .* normcdf(-(phi2-rho*phi1)./sqrt(1-rho^2)) ./ mvncdf([phi1,-phi2],[0,0],[1 -rho;-rho 1]) ...
    + mvnpdf([-phi1,phi2],[0,0],[1 -rho;-rho 1]) * normpdf(-phi1) .* normcdf((phi2-rho*phi1)./sqrt(1-rho^2)) ./ mvncdf([-phi1,phi2],[0,0],[1 -rho;-rho 1]) ...
    - mvnpdf([-phi1,-phi2],[0,0],[1 rho;rho 1]) * normpdf(-phi1) .* normcdf(-(phi2-rho*phi1)./sqrt(1-rho^2)) ./ mvncdf([-phi1,-phi2],[0,0],[1 rho;rho 1]);
end