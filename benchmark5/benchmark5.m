% Benchmark5 
% Magnetoquasistatic T-Omega formulation of inductor

load('benchmark5.mat');


% Setup matrices and rhs
Mmu = Tomega.Mmu;
Ss = Tomega.Ss;
Mrho = Tomega.Mrho;
C = Tomega.C;
idxdof_t = Tomega.idxdof_t;
np = Tomega.np_t;
Xs_t = Tomega.Xs_t;

Mt = [ 	sparse(np,np)		sparse(np, 3*np)	;	
			Mmu*Ss'				Mmu					;	
		];	
Kt = [	Ss*Mmu*Ss'			Ss*Mmu				;
			sparse(3*np, np)	C'*Mrho*C			;
		];
		
Mt_red = Mt(idxdof_t, idxdof_t);
Kt_red = Kt(idxdof_t, idxdof_t);

f=500;
omega = 2*pi*f;
Ampl = 1;
f_t = @(t)	[	-Ss*Mmu*Xs_t*Ampl*sin(omega*t);
					-Mmu*Xs_t*Ampl*omega*cos(omega*t); 
				](idxdof_t);


% Time stepping
T0 = 0;
Tend = 2e-3;
dt = 2e-5;
t = T0:dt:Tend;
n_t = length(t);

% Implicit Euler
xtomega = sparse(length(idxdof_t),n_t);
for i=1:n_t-1
	fprintf(['T-Omega time step t = ' num2str(t(i+1)) '\n']);
	xtomega(:,i+1) = (Mt_red/dt+Kt_red)\(f_t(t(i+1)) + Mt_red/dt*xtomega(:,i));
end %for

% Reconstruct solution
t_to = sparse(3*np, n_t);
psi_to = sparse(np, n_t);
psi_to(idxdof_t(1:np-1), :) = xtomega(1:np-1, :);
t_to(idxdof_t(np:end)-np, :) = xtomega(np:end, :);
h= Xs_t*sin(omega*t) + t_to + Ss'*psi_to;
	
Wmagn_to = 1/2*sum(h.*(Mmu*h),1);

% Plot 	
figure;
plot(t, Wmagn_to, 'b-');
xlabel('time (s)');
ylabel('Magnetic Energy (J)');

% Store results
dlmwrite('benchmark5.csv',{'t','W'})
dlmwrite('benchmark5.csv',[t(:) Wmagn_to(:)],'-append')