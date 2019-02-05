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
fprintf ('thermal velocity');
V=sqrt(C.kb*T/m)

%meanfreepath
fprintf ('mean free path');
MFP=V*Tmn

%assign location
fprintf ('assigned particles');

%choose how many Particles u want!!!!!!!!!!!:
p=5;


p
x=200e-9*rand(p,1);
y=100e-9*rand(p,1);
o=2*pi*rand(p,1);
Vx=cos(o)*V;
Vy=sin(o)*V;

dt=2e-14;%time interval
xdt=Vx*dt;
ydt=Vy*dt;
y1=y;
x1=x;


step=100;
n=1;
ss=1;
for t=1:1:step
    x=x+xdt;
    y =y+ydt;
    
    x1(x >l) = - (l - x1(x >l));
    x(x > l) = x(x > l)-(l);
    x1(x <0)    = l - x1(x <0);
    x(x < 0)    = x(x < 0 )+(l);
    
    ydt(y > w) = - ydt(y > w);
    y(y > w) = (w)-(y(y > w)-(w));
    ydt(y < 0) = -ydt(y < 0 );
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
    y1=y;
    x1=x;
    pause(0.1);
    hold on;
    
    %temperature plot
    figure (2)
    Vparticle = sqrt(Vx.^2+Vy.^2);
    TParticle = (m*Vparticle.^2)/(C.kb);
    ave=(sum(TParticle)/p);
    scatter(t,ave)
    xlim ([0, step]);
    ylim ([0, 650]);
    xlabel("steps")
    ylabel("K")
    title("Temperature plot")
    hold on
 
end