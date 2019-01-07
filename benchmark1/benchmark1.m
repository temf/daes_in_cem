% Benchmark3
% Spiral inductor

load('benchmark1.mat');

% Setup matrices and rhs
A = [sparse(3*np,3*np),    -Mmui*C; ...
     Mepsi*C',            sparse(3*np,3*np)];
g = @(t) [sparse(3*np,1);...
          -Mepsi*j_space*signal(t)];
u0 = sparse(6*np,1);
f  = @(t,u) A*u + g(t);

% Init time stepping
mu0    = 4*pi*1e-7; 
eps0   = 8.8541878176e-12;
dt     = sqrt(mu0*eps0)/sqrt((1/min(diff(x)))^2+(1/min(diff(y)))^2+(1/min(diff(z)))^2);
t      = 0:dt:2e-11;

% Init fields
h_new   = zeros(3*np,1); 
e_new   = zeros(3*np,1); 
h_old   = zeros(3*np,1); 
e_old   = zeros(3*np,1); 

% Set outputs
n     = 100;
t_out = zeros(1,ceil(length(t)/n));
h     = zeros(3*np,length(t_out)); 
e     = zeros(3*np,length(t_out));

% Output voltage
v = zeros(1,length(t));

% Leapfrog
for i=1:length(t)
    t_lf = t(i);
    js = j_space*signal(t_lf-dt/2);
    e_new = e_old+dt*(Mepsi*(C'*h_old-js));
    h_new = h_old-dt*(Mmui*(C*e_new));
    if mod(i-1,n)==0
        fprintf('Store step %d/%d at time: %e s\n',i,length(t),t_lf);
        k=(i-1)/n+1;
        e(:,k) = e_new;
        h(:,k) = h_new;
        t_out(k)= t_lf;
    end
    v(i) = -e_new(find(j_space));
    h_old = h_new;
    e_old = e_new;
end

% Plot output voltage
figure;
iv = signal(t);
plot(t,iv,t,v);
legend('Current','Voltage');
xlim([0 2e-11]);
ylabel('current (I) and voltage (V)');
xlabel('time (s)');
grid on

% Store results
dlmwrite('benchmark1.csv',{'t','v','i'})
dlmwrite('benchmark1.csv',[t(:) iv(:) v(:)],'-append')