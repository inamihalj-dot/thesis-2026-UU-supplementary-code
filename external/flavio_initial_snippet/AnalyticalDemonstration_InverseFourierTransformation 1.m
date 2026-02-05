clear all
close all
clc

% Frequencies
Nf=16;
fmin=250; % [MHz]
fmax=1300;% [MHz]
f=linspace(fmin,fmax,Nf)*10^6;
w=2*pi*f;

% Space Domain
x0=0.25; % Position of the perturbation [m]
dr = 0.001;
r=dr:dr:1;
r0=find(r==x0);

% Electrical Properties 
eps0 = 8.854187817 * 10^-12;
mu0 = 4*pi*10^-7;

eps_0 = eps0*56; % Muscle
%sigma_0 = 0.88;
sigma_0 = 0; 

eps_1 = eps0*11; % Fat
%sigma_1 = 0.1;
sigma_1 = 0; 


% Complex Permittivity 
epsilon0 = repmat(eps_0 + sigma_0./(1i.*w),length(r),1);
epsilon1 = epsilon0;
epsilon1(r0,:) = eps_1 + sigma_1./(1i.*w); % Perturbation

% Complex wavenumber 
k = w .* sqrt(mu0.*(eps_0+1i.*(sigma_0./w)));
beta = real(k);
alfa = imag(k);

% E-Field - Cylindrical Wave Approxximation
E=[];
for ii=1:Nf                                    
    E(:,ii) = 1./sqrt(r).*exp(1i.*k(ii).*r);
end


% Beaverstone Equation
E1E0 = E.*E;
Delta_Epsilon = epsilon1-epsilon0;
dS = Delta_Epsilon.*E1E0;
dS(isnan(dS))=0;
dTdz = repmat((2.*pi.*r(:)),1,Nf); % Integrating over the entire surface of the cylinder of radius r

% Delta S-Parameter
DS = 1i.*w.*sum(dS.*dr.*dTdz); 

% Discrete Inverse "Fourier" Transformation
dw=w(2)-w(1);
for kk=1:length(r)
    for ii=1:Nf       
        Psi(kk,ii) = (sqrt(mu0*eps_0)/pi)*(DS(ii)./(1i.*w(ii))).*exp(-1i.*2.*beta(ii).*r(kk)).*dw;
    end
end

f1=figure,plot(r,abs(sum(Psi,2)),'LineWidth',2)
hold on,plot([r(r0) r(r0)],[0 max(abs(sum(Psi,2)))*1.05],'--r','LineWidth',2)
title(['$\Delta\epsilon(r = ',num2str(r(r0)),')\neq 0$'],'FontSize',18, 'Interpreter','latex','Color','r')
xlim([0 0.3])
xlabel('$r [m]$','FontSize',16, 'Interpreter','latex')
ylabel('$\Psi_{ii}(r)$','FontSize',16, 'Interpreter','latex')
saveas(f1,'Figure_Psi.jpg')
