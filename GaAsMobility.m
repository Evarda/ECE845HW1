clc;
clear variables;
close all;

% Constants
hbar=6.582119*10^-16; %[eV*s]
hbarJ=1.0545718*10^-34; %[J*s]
kb=8.6173303*10^-5; %[eV/K]
kbJ=1.38064852*10^-23; %[J/K]
T=300; %[K]
ep0=8.854187817*10^-12; %[F/m]
e=1.6021766208*10^-19; %[C]
m0=9.10938356*10^-31; %[kg]

Ec=0;

% GaAs
effm = 0.067*m0; %[kg] for gamma, L, X respectively
rho = 5.36/1000*(100^3); %[kg/m^3]
vs = 5.24*10^5/100; %[m/s]
epr0 = 12.90;
eprInf = 10.92;

nE=50;
E=linspace(0.000001,2,nE);

% Acoustic Phonon Scattering
Dac = 7.01; %[eV] for gamma, L, X respectively

% Polar Optical Phonon Scattering
E0 = 0.03536; %[eV]
w0 = E0/hbar; % [1/s]
N0=(exp(E0/(kb*T))-1)^(-1);

% Ionized Impurity Scattering
dNI = 100;
NI = 1e23; % [1/m^3]
Z = 1;

GammaMAcoustic(1:nE)=0;
GammaMIonImp(1:nE,1:length(NI))=0;
GammaMPop(1:nE)=0;
GammaTot(1:nE,1:length(NI))=0;

k(1:nE)=0;
Energy=linspace(0.0001,2,nE);

for i=1:nE
    
    k(i) = sqrt(2*effm*E(i)/(hbar*hbarJ));

% Density of States
g3dAcoustic = sqrt(2)/(pi^2*hbar^3)*effm^(3/2)*sqrt(E(i)-Ec);

% Acoustic Phonon Scattering
GammaMAcoustic(i) = 2*pi/(hbarJ*hbar)^(1/2)*Dac^2*kb*T/(2*rho*vs^2)*g3dAcoustic;

% Polar Optical Phonon Scattering
GammaMPop(i)=(e^2*w0*(epr0/eprInf-1))/(4*pi*epr0*ep0*sqrt(hbarJ/hbar)*hbarJ*sqrt(2*E(i)/effm(1)))* ...
    (N0*sqrt(1+E0/E(i))+(N0+1)*sqrt(1-E0/E(i))-E0*N0/E(i)*asinh(sqrt(E(i)/E0))+E0*(N0+1)/E(i)*asinh(sqrt(E(i)/E0-1)));

GammaMPop=real(GammaMPop);
    
% Ionized Impurity Scattering
Ld=sqrt(ep0*eprInf*kbJ*T/(e^2*NI)); %[m]
gamma=sqrt(8*effm*E(i)*Ld^2/(hbar*hbarJ));
GammaMIonImp(i)=(hbar/hbarJ)^(3/2)*(NI*e^4)/(16*sqrt(2*effm(1))*pi*eprInf^2*ep0^2)*(log(1+gamma^2)-gamma^2/(1+gamma^2))*E(i)^(-3/2);

GammaTot(i) = GammaMAcoustic(i)+GammaMIonImp(i)+GammaMPop(i);

end
    
% Numerical Integration

ExpTauEAcoustic = 0;
ExpEAcoustic = 0;   

ExpTauEIonImp = 0;
ExpEIonImp = 0; 

ExpTauEPop = 0;
ExpEPop = 0; 

ExpTauETot = 0;
ExpETot = 0;

for i=1:nE
    f0=exp(-Energy(i)/(kb*T));
    g3d = sqrt(Energy(i));
    
    ExpTauEAcoustic = ExpTauEAcoustic + 1/GammaMAcoustic(i)*Energy(i)*f0*g3d;
    ExpEAcoustic    = ExpEAcoustic + Energy(i)*f0*g3d;
    
    ExpTauEIonImp = ExpTauEIonImp + 1/GammaMIonImp(i)*Energy(i)*f0*g3d;
    ExpEIonImp    = ExpEIonImp + Energy(i)*f0*g3d;
    
    ExpTauEPop = ExpTauEPop + 1/GammaMPop(i)*Energy(i)*f0*g3d;
    ExpEPop    = ExpEPop + Energy(i)*f0*g3d;
    
    ExpTauETot = ExpTauETot + 1/GammaTot(i)*Energy(i)*f0*g3d;
    ExpETot    = ExpETot + Energy(i)*f0*g3d;
    
end

% Find Tau
AvgTauAcoustic  = ExpTauEAcoustic/ExpEAcoustic;
AvgTauIonImp  = ExpTauEIonImp/ExpEIonImp;
AvgTauPop  = ExpTauEPop/ExpEPop;
AvgTauTot  = ExpTauETot/ExpETot;

% Calculate Mobilities
MobilityAcoustic = 100^2*e*AvgTauAcoustic/effm; % cm^2/(Vs)
MobilityIonImp = 100^2*e*AvgTauIonImp/effm; % cm^2/(Vs)
MobilityPop = 100^2*e*AvgTauPop/effm; % cm^2/(Vs)

% Calculate Total Mobility
MobilityTot = 100^2*e*AvgTauTot/effm; % cm^2/(Vs)

% Matthiessen's Rule
Mobility = (1/MobilityAcoustic + 1/MobilityIonImp + 1/MobilityPop).^(-1); % cm^2/(Vs)

disp(['The acoustic mobility is ' num2str(MobilityAcoustic)])
disp(['The ionized impurity mobility is ' num2str(MobilityIonImp)])
disp(['The POP mobility is ' num2str(MobilityPop)])
disp(['The mobility from Matthiessens Rule is ' num2str(Mobility)])
disp(['The mobility from total momentum relaxation rate is ' num2str(MobilityTot)])

% Conductivities
conductivity = NI*e*Mobility;
conductivityTot = NI*e*MobilityTot;

disp(['Conductivity is ' num2str(conductivity)])
disp(['Conductivity is ' num2str(conductivityTot)])