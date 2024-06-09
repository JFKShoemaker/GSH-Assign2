function [phi] = flexural (n,D_M3, rho_c, rho_m, g, R)

    phi = ( 1+ (D_M3/((rho_m-rho_c)*g))*((2*n+1)/(2*R))^4) ^ -1;

end
