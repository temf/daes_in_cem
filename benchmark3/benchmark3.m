% Benchmark3
% High-voltage cable termination 

load('benchmark3.mat');

% Setup matrices
K = Ksigma(idxdof,idxdof);
M = Keps(idxdof,idxdof);
f1 = -sum(Ksigma(idxdof,idxbcs),2);
f2 = -sum(Keps(idxdof,idxbcs),2);

% Compute steady state
Uapp = 1000;
y0   = pcg(K,f1*Uapp,1e-6,1000);

% Setup rhs
tswitch = 1e-3;
pp      = interp1([0 tswitch 2*tswitch 1],[1 1 -1 -1],'linear','pp');
ftU     = @(t) Uapp*ppval(pp,t);
ppd     = mkpp([0 tswitch 2*tswitch 1],[0;-1/(2*tswitch);0]);
ftdUdt  = @(t) Uapp*ppval(ppd,t);
rhs     = @(t) f1*ftU(t)+f2*ftdUdt(t);

% Init time stepping
dt      = 1e-5;
t       = 0:dt:5e-3;
x       = zeros(length(y0),length(t));
x(:,1)  = y0;

% Euler
[L, U, P, Q] = lu (M/dt+K);
for j=2:length(t)
	b=(M/dt)*x(:,j-1) + rhs(t(j));
	x(:,j) = Q*(U\(L\(P*b)));
end

% Reconstruct solution
Phi = zeros(size(Ksigma,1),length(t)); 
Phi(idxdof,:) = x; 
Phi(idxbcs,:) = ones(length(idxbcs),1)*ftU(t);

% Plot
figure; 
v = ftU(t);
plot(t,v,'-'); 
xlabel('time (s)'); 
ylabel('Excitation (V)');

figure; 
W = 1/2*sum(Phi.*(Keps*Phi),1);
plot(t,W,'-b'); 
xlabel('time (s)'); 
ylabel('Energy (W)');

figure;
q = Keps*Phi;
Q = sum(q(idxbcs,:),1);
plot(t,Q,'-r'); 
xlabel('time (s)'); 
ylabel('charge at the high-voltage electrode (C)');

% Store results
dlmwrite('benchmark3.csv',{'t','v','W','Q'})
dlmwrite('benchmark3.csv',[t(:) v(:) W(:) Q(:)],'-append')