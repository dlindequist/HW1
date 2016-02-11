function F = root_MP(theta)

global n P_T z_T k beta mu alph delta b A 

for i=1:n
F(i) = theta(i).^(alph)-k/(beta*A*sum(P_T(i,:)*((1-mu)*(z_T-b)-k.*mu.*theta-((1-delta)*k)./(A*theta.^(-alph)))));   
end
