% Benchmark4
% Magnetoquasistatic A* formulation of inductor

load('benchmark4.mat');

% Setup matrices and rhs
C = Astar.C;
SS = Astar.Ss;
Msigma = Astar.Msigma;
Mnu = Astar.Mnu;
Xstr = Astar.Xstr;	
idxdof_a = Astar.idxdof_a;
np		 = Astar.np_a;
	
Ma = Msigma;
Ka = C'*Mnu*C;

Ma_red = Ma(idxdof_a, idxdof_a);
Ka_red = Ka(idxdof_a, idxdof_a);

f=500;
omega = 2*pi*f;
Ampl = 1;
f_a = @(t)	[	Xstr*Ampl*sin(omega*t);
			](idxdof_a);

% Time stepping
T0 = 0;
Tend = 2e-3;
dt = 2e-5;
t = T0:dt:Tend;
n_t = length(t);

				
% Implicit Euler				
xastar = sparse(length(idxdof_a),n_t);
for i=1:n_t-1
	fprintf(['Astar time step t = ' num2str(t(i+1)) '\n']);
	xastar(:,i+1) = pcg((Ma_red/dt+Ka_red), (f_a(t(i+1)) + Ma_red/dt*xastar(:,i)), 1e-5, 100);
end %for
	
% Reconstruct solution
a = sparse(3*np, n_t);	
a(idxdof_a, :) = xastar;
b = C*a;
	
Wmagn_a = 1/2*sum(b.*(Mnu*b),1);

% Plot
figure;
plot(t, Wmagn_a, 'b-');
xlabel('time (s)');
ylabel('Magnetic Energy (J)');
	
% Store results
dlmwrite('benchmark4.csv',{'t','W'})
dlmwrite('benchmark4.csv',[t(:) Wmagn_a(:)],'-append')