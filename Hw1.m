clear all

%% Tauchen's Method
global n m p z k beta mu alph delta b A 

n=50; % number of grid points
sigma_epsilon=0.0034;
rho_z=0.9895;
sigma_z=sigma_epsilon/(sqrt(1-(rho_z)^2));

lambda=3;
z_upper=lambda*sigma_z;
z_lower=-lambda*sigma_z;

z=zeros(n,1); % creating equally distant z's
for i=1:n
z(i)=z_lower+(z_upper-z_lower)/(n-1)*(i-1);
end

m=zeros(n-1,1); % creating midpoints
for i=1:(n-1)
    m(i)=(z(i+1)+z(i))/2;
end

p=zeros(n,n); % create Markov matrix
for i=1:n
p(i,1)=normcdf((m(1)-rho_z*z(i))/sigma_epsilon);
p(i,n)=1-normcdf((m(n-1)-rho_z*z(i))/sigma_epsilon);
for j=2:n-1
p(i,j)=normcdf((m(j)-rho_z*z(i))/sigma_epsilon)-normcdf((m(j-1)-rho_z*z(i))/sigma_epsilon);
    end
end

%% solving for theta

k=0.6;
beta=0.999;
mu=0.72;
b=0.4;
delta=0.0081;
A=0.158;
alph=0.72;


fun=@root_MP;
theta0=zeros(n,1);
theta=fsolve(fun,theta0)

%% Simulating from the Markov Chain

z_start=z(1);
z_sim=zeros(n,1);
z_sim(1)=z_start;

location=1;
for i=2:n
   uniform_random=rand;
   for j=1:50
   if uniform_random<sum(p(location,j)), break;
   end
   location=find(uniform_random<sum(p(location,j)));
   end
   z_sim(i)=z(location);
   end
