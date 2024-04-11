% This Matlab code computes the ground state energy for helium (He) atom
% using the Car-Parrinello approach in terms of the quantum
% molecular dynamics method (qmd) [1,2]. 
%
% The core Hamiltonian matrix (H_core), overlap matrix (S_ov) and two-electron integrals (tei) (He_6_31G_tei.txt) are computed 
% by my own developing code, using the 6_31G basis set. An obtained total energy is compared with
% that from the restricted self-consistent field (SCF) calculation [2]. 
%
% Ref. [1]: Chapter 9 of "Computational Physics" book by Jos Thijssen. 
% Ref. [2]: https://github.com/aromanro/PythonCompphys/blob/master/Car-Parrinello.ipynb
% Ref. [3]: A. Szabo and N. S. Ostlund "Modern Quantum Chemistry" book 
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
%
% April 11, 2024 & University of North Dakota 
%
function [] = He_6_31g_qmd
clc;
%
format long
%
S_ov = [1.         0.63414774;
        0.63414774 1.        ];
%
H_core = [-1.87962624 -1.54640933;
          -1.54640933 -1.29518958];
%
%
dim = 2; %  (4s) -> [2s] = 2x1 = 2
%
itermax = 100; tol = 1e-6;
%
tei_n = 16;             % = 2^4, .i.e., all values of TEI
%
read_tei_data = fopen('He_6_31G_tei.txt', 'r');               % data of two-electron integral in atomic basis set
tei_data_n5 = textscan(read_tei_data, '%d %d %d %d %f');
%
p = zeros(tei_n,1); q = zeros(tei_n,1); r = zeros(tei_n,1); s = zeros(tei_n,1); vals = zeros(tei_n,1);
p(1:tei_n) = tei_data_n5{1};
q(1:tei_n) = tei_data_n5{2};
r(1:tei_n) = tei_data_n5{3};
s(1:tei_n) = tei_data_n5{4};
vals(1:tei_n) = tei_data_n5{5};
for i = 1:tei_n
    tei(p(i),q(i),r(i),s(i)) = vals(i);
%    tei(q(i),p(i),r(i),s(i)) = vals(i);    
%    tei(p(i),q(i),s(i),r(i)) = vals(i);    
%    tei(q(i),p(i),s(i),r(i)) = vals(i);   
    %
%    tei(r(i),s(i),p(i),q(i)) = vals(i);    
%    tei(s(i),r(i),p(i),q(i)) = vals(i);        
%    tei(r(i),s(i),q(i),p(i)) = vals(i);        
%    tei(s(i),r(i),q(i),p(i)) = vals(i);            
end
%
Q_tei = tei;
%
C = [0.5, 0.5]';
C_old = C;
%
dt = 0.1; % time step
dt2 = dt * dt;
dt4 = dt2 * dt2;
gamma = 1.0; %frictional constant
%
mass = 1.00;
massMinusGamma = mass - 0.5*gamma*dt;
massPlusGamma = mass + 0.5*gamma*dt;
%
En_0_old = 0.;
%
for iter = 1:itermax
    iter
    %
    F = H_core;
    for p = 1:dim
        for q = 1:dim
            for r = 1:dim
                for s = 1:dim
                    F(p,q) = F(p,q) + C(r)*C(s) * (2.*Q_tei(p,q,r,s) - Q_tei(p,r,s,q));
                end
    
            end
    
        end
    end
    % 
    En_0 = sum(C.* (H_core + F) * C);
    %
    if (abs(En_0_old - En_0) < tol)
        break
    end
    En_0_old = En_0;
    %
    Ct = (2. *mass * C - massMinusGamma * C_old - (4/1) * F * C * dt2 )./massPlusGamma;
    %
    a = sum(S_ov * S_ov * C.* S_ov * C) * dt4;
    b = -2.* sum(S_ov * C.* S_ov * Ct) * dt2;
    c = sum(S_ov * Ct.* Ct) - 1.;
    %
    delta = b * b - 4. * a * c;
    %
    if (delta < 0.)
        break  
    end
    %
    sdelta = sqrt(delta);
    lam1 = (-b - sdelta)/(2.*a) ;
    lam2 = (-b + sdelta)/(2.*a) ;
    %
    if  ( lam1 < 0.  ) 
        lam = lam2;
    else
        lam = lam1;
    end
    Ct = Ct - lam * dt2 * S_ov * C;
    %
    C_old = C;
    C = Ct;
    %

    [iter, iter * dt, En_0];

end
%%%
En_0    % -2.855159516975645au vs -2.855160422622062au from the SCF calculation



%%%
return
end
