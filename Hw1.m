clear all

global n m P_T P_R z_T z_R k beta mu alph delta b A 

n=50; % number of grid points
sigma_epsilon=0.0034;
rho_z=0.9895;
sigma_z=sigma_epsilon/(sqrt(1-(rho_z)^2));

%% Tauchen's Method

lambda=3;
z_upper=lambda*sigma_z;
z_lower=-lambda*sigma_z;

z_T=zeros(n,1); % creating equally distant z's (z-grid)
for i=1:n
z_T(i)=z_lower+(z_upper-z_lower)/(n-1)*(i-1);
end

m=zeros(n-1,1); % creating midpoints
for i=1:(n-1)
    m(i)=(z_T(i+1)+z_T(i))/2;
end

P_T=zeros(n,n); % create Markov matrix
for i=1:n
P_T(i,1)=normcdf((m(1)-rho_z*z_T(i))/sigma_epsilon);
P_T(i,n)=1-normcdf((m(n-1)-rho_z*z_T(i))/sigma_epsilon);
for j=2:n-1
P_T(i,j)=normcdf((m(j)-rho_z*z_T(i))/sigma_epsilon)-normcdf((m(j-1)-rho_z*z_T(i))/sigma_epsilon);
    end
end

%% Rouwenhorst's Method

p_R=(1+rho_z)/2;
q_R=p_R;
v=sqrt((n-1)/(1-rho_z^2))*sigma_z;

for i=1:n
    z_R(i)=-v+(2*v)/(n-1)*(i-1);    % creating the z-grid
end

P_R=[p_R 1-p_R;1-q_R q_R];      % creating transition matrix according to algorithm
for i=1:(n-2)
zeros_R=zeros(i+1,1);
P_R=p_R.*[P_R zeros_R;zeros_R' 0]+(1-p_R).*[zeros_R P_R; 0 zeros_R']+(1-q_R).*[zeros_R' 0; P_R zeros_R]+q_R.*[0 zeros_R'; zeros_R P_R];  
end

P_R=vertcat(P_R(1,:), P_R(2:49,:)/2, P_R(50,:));    % divide all but first and last row by 2


%% Solving for thetas

k=0.6;
beta=0.999;
mu=0.05;        % baseline: 0.72 / Hagedorn Manovski: 0.05

b=0.95;          % baseline: 0.4 / Hagedorn Manovski: 0.95
delta=0.0081;
A=0.158;
alph=0.72;


fun=@root_MP;
theta0=zeros(n,1);
theta=fsolve(fun,theta0);

%% Simulating from the Markov Chain

rounds=50;  % choose how many periods you wanna simulate       
start=5;      % pick gridpoint you wanna start at
z_start=z_T(start);
z_sim=zeros(n,1);
z_sim(1)=z_start;

rng('default');
rng(5);
uniform_random=rand(rounds,1);


location=zeros(rounds,1);
location(1)=start;
for i=2:rounds
   for j=1:50
   if uniform_random(i)>sum(P_T(location(i-1),1:j)) 
   else break 
   end
   end
   location(i)=j;
   z_sim(i)=z_T(location(i));
end

%% Obtaining endogenous variables from simulated productivity shocks

% get labor market tightness theta
theta_sim=zeros(rounds,1);
for i=1:rounds
theta_sim(i)=theta(location(i));
end

% get match probability of vacancy q
q_sim=zeros(rounds,1);
for i=1:rounds
q_sim(i)=A*(1/theta(location(i)))^(alph); 
end

% get match probability of worker p
p_sim=zeros(rounds,1);
for i=1:rounds
p_sim(i)=theta(location(i))*q_sim(location(i));    
end

% get unemployment rate unemp
unemp_sim=zeros(rounds,1);
for i=1:rounds
unemp_sim(i)=delta./(delta+theta_sim(location(i))*q_sim(location(i)));    
end

% plot endogenous variables
figure
subplot(2,2,1)
plot(z_sim)
title('simulated productivity shocks')

subplot(2,2,2)
plot(theta_sim)
title('simulated labor market tightness')

subplot(2,2,3)
h1 = plot(1:rounds, p_sim, 'blue');
hold on;
h2 = plot(1:rounds, q_sim, 'red');
legend([h1 h2],{'p_{sim}', 'q_{sim}'});
plot(1:rounds,p_sim,1:rounds, q_sim)
title('simulated matching probabilities')

subplot(2,2,4)
plot(unemp_sim)
title('simulated unemployment rate')
