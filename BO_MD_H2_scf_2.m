% This Matlab code computes the dynamics of internuclear separation (R) for
%  molecular hydrogen (H2) using the Born-Oppenheimer (BO) approach in terms of the quantum
% molecular dynamics method (qmd) [1, 2]. An internuclear separation (R)
% converges to R = 1.386 (au), which is close to exact value, 1.40 (an equilibrium separation).
% The self-consistent field (SCF) calculation is done using pseudo-spectral method by solving the Hartree-Fock equation [3]. 
%
% Solving equation: M * d^2/dt^2 = - nabla (Psi_0|H_el|Psi_0)
%                     H_el*Psi_0 = E0 * Psi_0    
% 
% Gradient of potential energy surface is computed by a finite difference scheme and
% velocity-Verlet method is used to integrate an equation [2]. 
%
% Ref. [1]: "Ab Initio Molecular Dynamics Basic Theory and Advanced Methods" by Dominik Marx and Jürg Hutter, Cambridge University Press (2009). 
% Ref. [2]: Chapter 9 of "Computational Physics" book by Jos Thijssen. 
% Ref. [3]: Tsogbayar Tsednee, PhD thesis at York University, 2014; 
%           https://yorkspace.library.yorku.ca/server/api/core/bitstreams/36f42c0d-c3bd-4a45-815b-812b1c736b60/content
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
%
% September 24, 2025 & Institute of Physics and Technology, Mongolian Academy of Sciences
%
function [] = BO_MD_H2_scf_2
clc;
% ----------------------------------------------------------------
%
gamma = 20.;
Mass_H = 1836.15;
dR = 0.01;
%
R0 = 1.3500;
dt = 5.0;
%
v0 = 0.;
%
N_step = 300;
%
fileID_save_data_1 = fopen('BO_MD_H2_scf_2.txt','w');
%
for ii = 1:N_step
    %
    [En_plus_dR] = H2_scf_plus_dR(R0, dR);
    [En_minus_dR] = H2_scf_plus_dR(R0, -dR);
    dE_R0 = (En_plus_dR - En_minus_dR)/(2*dR);
    RHS = dE_R0 - 1./R0^2; 
    a0 = -(1/Mass_H) * RHS ;
    v_half = v0 + 0.5 * dt * a0;
    R1 = R0 + (1./(1. + gamma*dt/(2*Mass_H))) * v_half * dt;    
    %
    % 2nd step
    [En_plus_dR1] = H2_scf_plus_dR(R1, dR);
    [En_minus_dR1] = H2_scf_plus_dR(R1, -dR);
    dE_R1 = (En_plus_dR1 - En_minus_dR1)/(2*dR);
    RHS_1 = dE_R1 - 1./R1^2; 
    a1 = -(1/Mass_H) * RHS_1 ;
    v1 = ((1. - gamma*dt/(2*Mass_H))/(1. + gamma*dt/(2*Mass_H))) * v_half + 0.5 * a1 * dt;    
    %
    v0 = v1;
    R0 = R1;    
    %
    KE = 0.25*Mass_H*v1^2;
    %
    [ii * dt, R1, v1, KE]
    output = [ii * dt, R1, v1, KE];    
    %
    fprintf(fileID_save_data_1, '%4.4f \t %4.12f \t %4.12f \t %8.12f\n', output); 


end
%
fclose(fileID_save_data_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
read_md_data = fopen('BO_MD_H2_scf_2.txt', 'r');               % 
read_md_data = textscan(read_md_data, '%f %f %f %f');
md_step_ii = read_md_data{1};
md_R = read_md_data{2};
md_vel = read_md_data{3};
md_KE = read_md_data{4};
%

figure(1)
hold on
plot(md_step_ii, md_R, 'b', LineWidth=1.2)
hold off
xlabel('\mbox{Time}\,(au)','Interpreter','latex') % ,'fontsize',16
ylabel('$R\,(au)$','Interpreter','latex') % , 'Rotation',0 ,'Rotation',1
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on

figure(2)
hold on
plot(md_step_ii, md_vel, 'b', LineWidth=1.2)
hold off
xlabel('\mbox{Time}\,(au)','Interpreter','latex') % ,'fontsize',16
ylabel('$Velocity$','Interpreter','latex') % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on


figure(3)
hold on
plot(md_step_ii, md_KE, 'b', LineWidth=1.2)
hold off
xlabel('\mbox{Time}\,(au)','Interpreter','latex') % ,'fontsize',16
ylabel('$Kinetic\, energy$','Interpreter','latex') % , 'Rotation',0
%axis([0. 3. 0.000 2.00])
set(gca,'FontSize',16)
box on


%%%
return
end


%%%
function [En_plus_dR] = H2_scf_plus_dR(R_int, dR)
%%%
Rint = R_int + dR; 
%
N = 16.; M = 8; a = 1.; b = 20.;    % Initial computaional parameters; you may change them
%Rint = 1.350000251866610 + 0.01; % internuclear separation
itermax = 100; tol = 10^(-12);  % max number of iteration and tolerance; you may change them
% mu part begins
[mu,wmu,D] = legDC2(N,a,b); mup = mu;
D2 = D*D;
Tmu = diag(mu.^2 - 1.)*(2/(b-a))^2*D2 + ...
      2.*diag(mu)*(2/(b-a))^1*D;
% mu part ends  
%%%
% nu part begins
[nu,wnu,Dn] = legDC2(M,-a,a); nup = nu;
D2n = Dn*Dn;
Tnu = diag(1. - nu.^2)*(2/(a-(-a)))^2*D2n - ...
      2.*diag(nu)*(2/(a-(-a)))*Dn;
% nu part ends

% right hand side of an equation

%%%
mup(1);
[mump,nunp] = meshgrid(mup,nup); 
mump = mump(:);
nunp = nunp(:);
Smunup = diag(mump.^2 - nunp.^2);
%%%

%%% initial rho0 begins ---
% LCAO approximation
sm_ov = (1.+Rint+(1./3.)*Rint^2)*exp(-Rint);
Nconst = 1./(sqrt(pi)*sqrt(2.+2.*sm_ov)); 
psi0 = Nconst*(exp(-0.5*Rint*(mump+nunp)) + ...  %%% +
               exp(-0.5*Rint*(mump-nunp)));
                            
rho0 = psi0.^2;
%rho0';
%%% initial rho0 ends ---

% Impose boundary conditions by replacing appropraite rows:
bmul = find(abs(mump) == mup(1) );
%
Hmn_p = (4./Rint^2)*(kron(Tmu(1:N+1,1:N+1),eye(M+1)) + ...
                     kron(eye(N+1),Tnu(1:M+1,1:M+1))); 

Hmn_p(bmul,:) = zeros(size(bmul,1),(N+1)*(M+1));
Hmn_p(bmul,bmul) = eye(size(bmul,1));

%%%  BC implementation begins
rhs = zeros((N+1)*(M+1),1);
rhs(bmul) = (mump(bmul) == mup(1)).*(1./((Rint/2.).*sqrt(mump(bmul).^2 +...
                                                         nunp(bmul).^2 -...
                                                         1. ) ));
rhsrsh = reshape(rhs,M+1,N+1);
%%%  BC implementation ends

%%% Poisson solver begins    
%%% right hand side calculation begins
f = -Smunup*4*pi*rho0;
frsh = reshape(f,M+1,N+1);
%%% right hand side calculation ends
ff = frsh + 1.*rhsrsh;
f = reshape(ff,(N+1)*(M+1),1);
u = Hmn_p\f;  
uu = reshape(u,M+1,N+1);
uup = uu';
%%% Poisson solver ends

%%% Electron-electron interaction calculation begins
%uv = reshape(u,Nn+1,N+1)
%%% initial rho0 begins
psi0p = Nconst*(exp(-0.5*Rint*(mump+nunp)) + ... %%% +
                exp(-0.5*Rint*(mump-nunp)));
                            
rho0p = psi0p.^2;
%%% initial rho0 ends
rho0rshp = reshape(rho0p,M+1,N+1)';
%%%
Vee = 0.;
for i = 1:N+1
    for j = 1:M+1
        Vee = Vee + wmu(i)*wnu(j)*rho0rshp(i,j)*uup(i,j)*...
                    (Rint^3/8)*(mup(i)^2-nu(j)^2)*2*pi;
    end
end
Vee;
%%% Electron-electron interaction calculation ends

%%%
mu = mu(2:N+1);
[mum,nun] = meshgrid(mu,nu);
mum = mum(:);
nun = nun(:);
Smunu= diag(mum.^2 - nun.^2);

% right hand side of equation

%%% H2+ matrix element calculation begins
Tmunu = (4./Rint^2)*(kron(Tmu(2:N+1,2:N+1),eye(M+1)) + ...
                     kron(eye(N+0),Tnu(1:M+1,1:M+1)));             
Vcmunu = -(4./Rint)*diag(mum);
Hmn2 = -0.5*Tmunu + Vcmunu;
%Smunu= diag(mum.^2 - nun.^2);
%%% H2+ matrix element calculation ends
[Veci,Eni] = eig(Hmn2,Smunu);                                     % Eigenvalue problem
Eni = diag(Eni);
[foo, ij] = sort(Eni);
Eni = Eni(ij);
%[Eni(1),Eni(2),Eni(3),Eni(4),Eni(5)]  
%%%

%%% H2 energy calculation begins --------------------- ---------------------

%%% H2+ matrix element calculation begins
Tmunuh2 = (4./Rint^2)*(kron(Tmu(1:N+1,1:N+1),eye(M+1)) + ...
                     kron(eye(N+1),Tnu(1:M+1,1:M+1)));             
Vcmunuh2 = -(4./Rint)*diag(mump);
Hmnh2 = -0.5*Tmunuh2 + Vcmunuh2;
%Smunu= diag(mum.^2 - nun.^2);
%%% H2+ matrix element calculation ends

Hmnh2(bmul,:) = zeros(size(bmul,1),(N+1)*(M+1));
Hmnh2(bmul,bmul) = eye(size(bmul,1));

%%%  BC implementation begins
rhswf = zeros((N+1)*(M+1),1);
%rhswf(bmul) = (mump(bmul) == mup(1)).*(2.*1./(Rint.*mump(bmul)));
rhswf(bmul) = (mump(bmul) == mup(1)).*exp(-sqrt(-2.*Eni(1))*...    %%%Eni(1)
                                          (Rint/2.).*mump(bmul) );
%rhsrshwf = reshape(rhswf,Nn+1,N+1);
rhsrshwf = diag(rhswf);
%%%  BC implementation ends

%%% right hand side calculation begins
fh2 = Smunup;
Smunuh2 = fh2 - rhsrshwf ;
%%% right hand side calculation ends

%%% Hartree-Fock hamiltonian: H2 energy calculation begins
Hmn = 1*Hmnh2 + diag(u)*Smunup;
[Vec,En] = eig(Hmn,Smunuh2);                                     % Eigenvalue problem
En = diag(En);
[foo, ij] = sort(En);
En = En(ij);
%[En(1),En(2),En(3),En(4),En(5)] ;

En_tot = 2*En(1) - Vee + 1/Rint ;
%%% H2 energy calculation ends --------------------- ---------------------

%%% wave function calculation begins
Vec = Vec(:,ij);  % The unnormalized eigenfunctions
V1 = Vec(:,1);    % The unnormalized eigenfunction for the ground state
%
V1 = reshape(V1,M+1,N+1)';
sm = 0.;
for i = 1:N+1
    for j = 1:M+1
        sm = sm + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                  (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;        
    end
end
sm; 
n_c = 1./sqrt(sm);      % The normalization constant
V1 = n_c*V1;
%figure(2)
%mesh(mup,nu,V1')
%%%
n_wf = 0.;
for i = 1:N+1
    for j = 1:M+1
%        V1(i,j) = n_c*V1(i,j);
        n_wf = n_wf + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                  (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;
    end
end
n_wf;  % check normalization 
rho0d = V1.^2;
%%% wave function calculation ends

%Iflag = 0.;  % converfence flag, 1 indicates convergence

for iter = 2:itermax
    
    Vee_old = Vee;

    iter; % 2
%    
    rho0p = rho0d;
    rho0pn = rho0p(1:N+1,1:M+1)';
    rho0pn = reshape(rho0pn,(N+1)*(M+1),1);
    rho0 = rho0pn;

%%% Poisson solver begins    
%%% right hand side calculation begins
    f = -Smunup*4*pi*rho0;
    frsh = reshape(f,M+1,N+1);
%%% right hand side calculation ends
%frsh(1:Nn+1,1) = 0.;
    ff = frsh + rhsrsh;
    f = reshape(ff,(N+1)*(M+1),1);
    %figure(1), clf, spy(Hmn_p), drawnow
    u = Hmn_p\f;   
    uu = reshape(u,M+1,N+1);
    uup = uu';
%%% Poisson solver ends

%%% Electron-electron interaction calculation begins
    rho0rshp = rho0d;
%%%
    Vee = 0.;
    for i = 1:N+1
        for j = 1:M+1
            Vee = Vee + wmu(i)*wnu(j)*rho0rshp(i,j)*uup(i,j)*...
                        (Rint^3/8)*(mup(i)^2-nu(j)^2)*2*pi;
        end
    end
%    Vee;
%%% Electron-electron interaction calculation ends

%%%
%%% H2 energy calculation begins

%%% Hartree-Fock hamiltonian: H2 energy calculation begins
    Hmn = Hmnh2 + diag(u)*Smunup;
    [Vec,En] = eig(Hmn,Smunuh2);                                     % Eigenvalue problem
    En = diag(En);
    [foo, ij] = sort(En);
    En = En(ij);
%    [En(1),En(2),En(3),En(4),En(5)] ; % orbital energy for ground state: En(1)
    En_tot = 2*En(1) - Vee + 1/Rint;
%%% H2 energy calculation ends

%%% wave function calculation begins
    Vec = Vec(:,ij);  % The unnormalized eigenfunctions
    V1 = Vec(:,1);    % The unnormalized eigenfunction for the ground state, 

    %reshape(V1,Nn+1,N)         ; 
    V1 = reshape(V1,M+1,N+1)';
    sm = 0.;
    for i = 1:N+1
        for j = 1:M+1
%        sm = sm + V1(i,j)*V1(i,j);
            sm = sm + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                      (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;        
        end
    end
    n_c = 1./sqrt(sm);      % The normalization constant
    V1 = n_c*V1;
%
    rho0ds = conj(V1).*V1;
%%% wave function calculation ends

    if (abs(Vee - Vee_old) < tol)
        break
    end
%   
    rho0d = rho0ds;
%   
end
%%%

%%%% some simple calculation for the ground state ---
% n_wf: check normalization of wave function
% smrsq: calculation of < |r^2| >  at given R to check wave function
% smQ: calculation of < |Q| >  at given R to check wave function 
n_wf = 0.;  
smrsq = 0.;
smQ = 0.;
for i = 1:N+1
    for j = 1:M+1
%        V1(i,j) = n_c*V1(i,j);
        n_wf = n_wf + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                  (Rint^3/8)*(mup(i).^2 -nu(j).^2)*2*pi;
%              
        smrsq = smrsq + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                    (Rint^2/4)*(mup(i).^2 + nu(j).^2 - 1)*...
                    (Rint^3/8)*(mup(i).^2-nu(j).^2)*2*pi;              

%%% <Q>---        
        smQ = smQ + wmu(i)*wnu(j)*V1(i,j)*V1(i,j)*...
                        0.5*(3.*(Rint^2/4)*mup(i).^2*nu(j).^2 - ...
                            (Rint^2/4)*(mup(i).^2 + nu(j).^2 - 1))*...
                            (Rint^3/8)*(mup(i).^2-nu(j).^2)*2*pi;
%%% <Q>---                  
                
    end
end
%n_wf;
%smrsq;   % < |r^2| >  at given R to check wave function 
%smQ;     % < |Q| >  at given R to check wave function 

% Output: 
[Rint, En(1), Vee, 2*En(1) - Vee, En_tot, smrsq, smQ]; % 
% 1.400000000000000  -0.594658571104492   0.658598144960419  -1.133629572883688   2.573929849836470   0.243288986798109 % 

En_plus_dR = 2*En(1) - Vee;
%%%

return
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
xi=(a*(1-x)+b*(1+x))/2;      % added by Tsogbayar Tsednee

% Compute the weights
w=(b-a)./(N*N1*P(:,N1).^2);  % added by Tsogbayar Tsednee
return
end
