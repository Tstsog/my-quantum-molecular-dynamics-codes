% This Matlab code computes the dynamics of internuclear separation (R) for
% hydrogen molecular ion (H2+) using the Born-Oppenheimer (BO) approach in terms of the quantum
% molecular dynamics method (qmd) [1, 2]. An internuclear separation (R)
% converges to R = 1.99972 (au), which is an equilibrium separation. 
%
% Solving equation: M * d^2/dt^2 = - nabla (Psi_0|H_el|Psi_0)
%                     H_el*Psi_0 = E0 * Psi_0    
% 
% Gradient of potential energy surface is computed by a finite difference scheme and
% velocity-Verlet method is used to integrate an equation [2]. 
%
% Ref. [1]: "Ab Initio Molecular Dynamics Basic Theory and Advanced Methods" by Dominik Marx and Jürg Hutter, Cambridge University Press (2009). 
% Ref. [2]: Chapter 9 of "Computational Physics" book by Jos Thijssen. 
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
%
% April 20, 2024 & University of North Dakota 
%
function [] = BO_MD_H2plus_ion_2
%
clc;
%
gamma = 15.;  % frictional (damping) constant 
Mass_H = 1836.15/2; % reduced mass (mu = mass_of_proton/2)
dR = 0.01;
%
R0 = 1.000;
dt = 10.0;
%
v0 = 0.;
%
N_step = 200;
%
fileID_save_data_1 = fopen('BO_MD_H2plus_ion_2.txt','w');
%
for ii = 1:N_step
    %
    [En_plus_dR] = H2plus_eig_values_for_sigma_states(R0, dR);
    [En_minus_dR] = H2plus_eig_values_for_sigma_states(R0, -dR);
    dE_R0 = (En_plus_dR - En_minus_dR)/(2*dR);
    RHS = dE_R0 - 1./R0^2; 
    a0 = -(1/Mass_H) * RHS ;
    v_half = v0 + 0.5 * a0 * dt - (gamma/Mass_H) * v0 * dt;
    R1 = R0 + v_half * dt;
    %
    % 2nd step
    [En_plus_dR1] = H2plus_eig_values_for_sigma_states(R1, dR);
    [En_minus_dR1] = H2plus_eig_values_for_sigma_states(R1, -dR);
    dE_R1 = (En_plus_dR1 - En_minus_dR1)/(2*dR);
    RHS_1 = dE_R1 - 1./R1^2; 
    a1 = -(1/Mass_H) * RHS_1 ;
    v1 = v_half + 0.5 * a1 * dt;
    %
    v0 = v1;
    R0 = R1;    
    %
    KE = 0.5*Mass_H*v1^2;
    %
    [ii * dt, R1, v1, KE]
    output = [ii * dt, R1, v1, KE];    
    %
    fprintf(fileID_save_data_1, '%4.4f \t %4.12f \t %4.12f \t %8.12f\n', output); 

end
%
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_md_data = fopen('BO_MD_H2plus_ion_2.txt', 'r');               % 
read_md_data = textscan(read_md_data, '%f %f %f %f ');
md_step_ii = read_md_data{1};
md_R = read_md_data{2};
md_vel = read_md_data{3};
md_KE = read_md_data{4};
%

figure(1)
plot(md_step_ii, md_R, 'b', LineWidth=1.2)
xlabel('\mbox{Time}\,(au)','Interpreter','latex') % ,'fontsize',16
ylabel('$R\,(au)$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on

%
figure(2)
hold on
plot(md_step_ii, md_vel, 'b', LineWidth=1.2)
hold off
xlabel('\mbox{Time}\,(au)','Interpreter','latex') % ,'fontsize',16
ylabel('$Velocity$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on

figure(3)
hold on
plot(md_step_ii, md_KE, 'b', LineWidth=1.2)
hold off
xlabel('\mbox{Time}\,(au)','Interpreter','latex') % ,'fontsize',16
ylabel('$Kinetic\, energy$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on


%%%
return
end


function [E0] = H2plus_eig_values_for_sigma_states(R0, dR)
%-----------------------------------------------------------------
% Compute eigenvalue for the sigma states for H2+ ion in prolate
% spheroidal coordinate
% Written by Tsogbayar Tsednee (PhD), Texas A&M University; 
% Contact: tsog215@gmail.com
% Feb 14, 2017
% ----------------------------------------------------------------
%clear;
%clc;
N = 24; M = 10.; a = 1.; b = 20.; % Initial computaional parameters; you may change them
%R = 2.0; % internuclear separation; you may change it
R = R0 + dR;
[mu,wrm,Dmu]=legDC2(N,a,b);  % Legendre diff. matrix and coordinate mu and weight
[nu,wrn,Dnu]=legDC2(M,-a,a); % Legendre diff. matrix and coordinate nu and weight
Dmu2 = Dmu*Dmu; Dnu2 = Dnu*Dnu; % kinetic energy matrix
Tmu = diag(mu.^2-1.)*(2/(b-a))^2*Dmu2 + 2.*diag(mu)*(2/(b-a))*Dmu;
Tnu = diag(1.-nu.^2)*Dnu2 - 2.*diag(nu)*Dnu;  
[mum,nun] = meshgrid(mu(2:N+1),nu); mum = mum(:); nun = nun(:);
Smunu= diag(mum.^2 - nun.^2);
Tmunu = (4./R^2)*(kron(Tmu(2:N+1,2:N+1),eye(M+1)) + ...
                     kron(eye(N),Tnu(1:M+1,1:M+1)));
Vcmunu = -(4./R)*diag(mum); % Coulomb potential 
Hmunu = -0.5*Tmunu + Vcmunu; % H0 hamiltonian matrix
  
En = sort(eig(Hmunu,Smunu));  % solve eigenvalues 
%[En(1),En(2),En(3),En(4),En(5)]; % eigenenergies for first 5 S-states
% Results
% -1.1024   -0.6676   -0.3609   -0.2554   -0.2358

E0 = En(1);

%%%
return
end

function [xi,w,D]=legDC2(N,a,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% legDc.m
%
% Computes the Legendre differentiation matrix with collocation at the 
% Legendre-Gauss-Lobatto nodes.
%
% Reference: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 05/26/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Truncation + 1
N1=N+1;

% CGL nodes
xc=cos(pi*(0:N)/N)';

% Uniform nodes
xu=linspace(-1,1,N1)';

% Make a close first guess to reduce iterations
if N<3
    x=xc;
else
    x=xc+sin(pi*xu)./(4*N);
end

% The Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;
while max(abs(x-xold))>eps

    xold=x;
        
    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
     
    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );
end

X=repmat(x,1,N1);
Xdiff=X-X'+eye(N1);

L=repmat(P(:,N1),1,N1);
L(1:(N1+1):N1*N1)=1;
D=(L./(Xdiff.*L'));
D(1:(N1+1):N1*N1)=0;
D(1)=(N1*N)/4;
D(N1*N1)=-(N1*N)/4;

% Linear map from[-1,1] to [a,b]
xi=(a*(1-x)+b*(1+x))/2;      

% Compute the weights
w=(b-a)./(N*N1*P(:,N1).^2);


%%%
return
end
