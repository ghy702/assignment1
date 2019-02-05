clc
clear
global C
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s²
C.am = 1.66053892e-27;
T=300;
l=200e-9;
w=100e-9;
m=0.26*C.m_0;
Tmn=0.2e-12;

%thermal velocity
V=sqrt(C.kb*T/m);


fprintf ('assigned particles');

%choose how many Particles u want!!!!!!!!!!!:
p=30;
p

RVx=randn(p,1)*sqrt(C.kb*T/m);
RVy=randn(p,1)*sqrt(C.kb*T/m);
RV=sqrt(RVx.^2+RVy.^2);

%Histogram
figure (3)
plot(hist(RV,10))
ylabel("scalar")
xlabel("bins")
title ("Histogram")

x=200e-9.*rand(p,1);
y=100e-9.*rand(p,1);
dt=2e-14;%time interval
xdt=RVx*dt;
ydt=RVy*dt;
y1=y;
x1=x;


MFPx = x;
MFPy = y;
MFP  = zeros(p,1);
MTBC = zeros(p,1);
scatter_number = zeros(p,1);
Pscat = 1-exp(-dt/Tmn); 


step=100;

for t=1:1:step
    PRand=rand(p,1);
    MFP(PRand<Pscat) = MFP(PRand<Pscat) + sqrt((xdt(PRand<Pscat)-MFPx(PRand<Pscat)).^2+(ydt(PRand<Pscat)-MFPy(PRand<Pscat)).^2);
    MTBC(PRand<Pscat) =  MTBC(PRand<Pscat) + sqrt((xdt(PRand<Pscat)-MFPx(PRand<Pscat)).^2+(ydt(PRand<Pscat)-MFPy(PRand<Pscat)).^2)./RV(PRand<Pscat);
    scatter_number(PRand<Pscat) = scatter_number(PRand<Pscat) + 1;
    MFPx(PRand<Pscat) = x(PRand<Pscat);
    MFPy(PRand<Pscat) = y(PRand<Pscat);
    RVxn = randn(p,1)*sqrt(C.kb*T/m);
    RVyn = randn(p,1)*sqrt(C.kb*T/m); 
    RVx(PRand<Pscat) = RVxn(PRand<Pscat); 
    RVy(PRand<Pscat) = RVyn(PRand<Pscat);
    RV=sqrt(RVx.^2+RVy.^2);
    
    xdt=RVx*dt;
    ydt=RVy*dt; 
    x=x+xdt;
    y =y+ydt;
    
    x1(x >l) = - (l - x1(x >l));
    x(x > l) = x(x > l)-(l);
    x1(x <0)    = l - x1(x <0);
    x(x < 0)    = x(x < 0 )+(l);
    
    RVy(y > w) = - RVy(y > w);
    y(y > w) = (w)-(y(y > w)-(w));
    RVy(y < 0) = -RVy(y < 0 );
    y(y < 0) = -y(y < 0);
    

    
    
    %2-D plot of particle trajectories
    figure (1)
    for n=1:1:p
        
    plot([x1(n), x(n)], [y1(n), y(n)],'r');

    end

    xlabel("length(200nm)")
    ylabel("width(100nm)")
    title ("2-D plot of particle trajectories")
    xlim([0 l]);
    ylim([0 w]);
    pause(0.1);
    hold on;
    y1=y;
    x1=x;
    
    %temperature plot
    figure (2)
    Vparticle = sqrt(RVx.^2+RVy.^2);
    TParticle = (0.5*m*Vparticle.^2)/(C.kb);
    ave(t)=(sum(TParticle)/p);
    plot(ave,'g')
    xlim ([0, step]);
    ylim ([0, 1000]);
    xlabel("steps")
    ylabel("K")
    title("Temperature plot")
    hold on
 
end
fprintf ('ave mean free path');
MFP = sum(MFP./scatter_number)/p
fprintf ('ave Tmn');
Tmnn = sum(MTBC./scatter_number)/p
