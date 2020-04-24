%PTC last updated 4/22/2020

%Example use of function "cylinderscatteringcalculatorfromopdf.m"
%Simple shear flow, Per=100/17
%Cylinder with R=30nm, L=500nm


clear all

%parameters for form factor
R=30;%cylinder radius (in Angstroms)
L=500;%cylinder length (in Angstroms)
volfrac=0.1;%volume fraction of cylinders
dsld=(6.33E-6)-(3.13E-6);%scattering length density difference between cylinder and solvent (in Angstroms^2)
bkgd=0.0001;%background scattering in cm^-1

% %Qx Qy bounds if not using lists, comment out if using manual Qx,Qy coordinates
Qxmin=-0.05; %qx lower bound in Angstrom^-1
Qxmax=0.05; %qx upper bound in Angstrom^-1
Qxnum=50; %number of pixels in x direction
Qxstep=(Qxmax-Qxmin)/(Qxnum-1);

Qymin=-0.05; %qy lower bound in Angstrom^-1
Qymax=0.05; %qy upper `bound in Angstrom^-1
Qynum=50; %number of pixels in y direction
Qystep=(Qxmax-Qxmin)/(Qxnum-1);

%initialize wavevectors where you want to make the calculations
nQ=Qxnum*Qynum;
Qx=zeros(nQ,1);
Qy=zeros(nQ,1);
Qz=zeros(nQ,1);
xcount=0;
for Qxi=Qxmin:Qxstep:Qxmax
    xcount=xcount+1;
    ycount=0;
    for Qyi=Qymin:Qystep:Qymax
        ycount=ycount+1;
        Qx(xcount+(Qxnum)*(ycount-1))=Qxi;
        Qy(xcount+(Qxnum)*(ycount-1))=Qyi;
    end
end


%flow and discretization parameters for numerical simulation of OPDF
nangle=20;%angle discretization
Dr=17; %same units as deformation rate
G=100;%Deformation Rate
Per=G/Dr;
AR=(L./(2*R))*sqrt(8*pi/(16.35*log(L/(2*R))));%Particle Effective Aspect Ratio
ARG=(AR.^2-1)/(AR.^2+1);%Shape Factor
tstart=0;
store=10000;%store results per this many timesteps
strain=max(100,0.1*pi*AR);%strain units
tstop=strain./G; %simulation time
%tstop=max(1/(Dr*AR^2),pi.*AR./G); %for large Per
dt=1000;%Only used for no flow calculations, otherwise keep large
%dt=0.01;%uncomment when in no flow

%Set initial condition
Yinitcond=ones(nangle)./(4*pi); %uniform (equilibrium) initial condition

%Velocity Gradient Tensor

% k = [k11 k12 k13]
%     [k21 k22 k23]
%     [k31 k32 k33]

k11=0;
k12=G;
k13=0;
k21=0;
k22=0;
k23=0;
k31=0;
k32=0;
k33=0;

k=[k11 k12 k13; k21 k22 k23; k31 k32 k33];

[Y,a2,a4,Ytot,phivect,thetavect] = FokkerPlanckDiluteRodSolver(nangle,k,Dr,ARG,Yinitcond,tstart,store,tstop,dt);

[I,c] = cylinderscatteringcalculatorfromopdf(Qx,Qy,Qz,Y,nangle,R,L,volfrac,dsld,bkgd);

%plotting
%simulated spectra
Ilog=log10(I);
figure
plot=scatter(Qx,Qy,4*150000/(nQ),Ilog,'sq','filled');
colormap jet
colorbar
h = colorbar;
ylabel(h, 'log(I(q))')
% xlim([Qxmin-Qxstep Qxmax+Qxstep])
% ylim([Qymin-Qystep Qymax+Qystep])
xlabel('$q_x [\AA^{-1}]$','Interpreter','latex','Fontname','TimesNewRoman')
ylabel('$q_y [\AA^{-1}]$','Interpreter','latex','Fontname','TimesNewRoman')
axis equal
set(gca,'FontSize',16)




