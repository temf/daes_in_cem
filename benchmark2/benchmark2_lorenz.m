% benchmark2: Lorenz gauge
% Metal bar

load('benchmark2_lorenz.mat')

% Setup matrices and rhs
n   = size(M,1);
x0  = zeros(n,1);
f1  = @(t) 2*pi*cos(2*pi*t);
f2  = @(t) sin(2*pi*t);
rhs = @(t) r*[f1(t); 0; f2(t); 0];

% Init time stepping
dt     = 1e-4;
t      = 0:dt:1;
x      = zeros(length(t),n);
x(1,:) = x0;

% Rescale equations
A = M/dt+K;
R = spdiags(1./max(abs(A'))',0,n,n);

% Euler
[L, U, P, Q] = lu (R*A);
for j=2:length(t)
	b=(M/dt)*x(j-1,:)' - rhs(t(j));
	x(j,:) = Q*(U\(L\(P*R*b)));
end

% Plot
figure;
v = f2(t);
plot(t,v,'-b');
xlabel('time (s)'); 
ylabel('voltage at the contact (V)');

figure;
iv = x(:,end-1);
plot(t,iv,'-r');
xlabel('time (s)'); 
ylabel('current through the contact (A)');

% Store results
dlmwrite('benchmark2_lorenz.csv',{'t','v','i'})
dlmwrite('benchmark2_lorenz.csv',[t(:) v(:) iv(:)],'-append')