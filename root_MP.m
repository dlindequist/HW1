function F = root_MP(theta)

global n p z k beta mu alph delta b A 

for i=1:n
F(i) = theta(i).^alph-k/(beta*A*sum(p(i,:)*((1-mu)*(z-b)-k.*mu.*theta-((1-delta)*k)./(A*theta.^(-alph)))));   
end

