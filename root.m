function F=root(theta)

for i=1:n
summing=sum(p(i,:)*((1-mu)*(z-b)-k*mu*theta-((1-delta)*k)./(A*theta.^(-alpha))));
F(i)=theta(i)^alpha-k/(beta*A*summing);      
end
