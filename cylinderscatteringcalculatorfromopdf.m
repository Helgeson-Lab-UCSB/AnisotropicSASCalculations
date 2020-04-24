function [I,c] = cylinderscatteringcalculatorfromopdf(Qx,Qy,Qz,Y,nangle,R,L,volfrac,dsld,bkgd)
%%Program solves for scattering from dilute cylinders oriented in a flow
%Run this after generating orientation distribution (Y)
%INPUTS
%Qx,Qy,Qz: arrays of the Qx,Qy,Qz coordinates where the scattering 
%    should be calculated in Angstroms^-1
%Y: The OPDF as a matrix generated from the FP simulations
%nangle: number of angles used as input to FP simulations
%R: radius of cylinder in inverse units of Q (Angstroms)
%L: length of cylinder in inverse units of Q (Angstroms)
%volfraction: volume fraction of cylinders
%dsld: scattering length density difference between cylinder and solvent
%bkgd: incoherent background scattering (in cm^-1)
%
%OUTPUTS
%I: Scattering intensity for Qx,Qy,Qz coordinates in cm^-1
%c: scaling factor for scattering intensity in cm^-1, subtract background
%   and divide by this to get the form factor scattering, 
%   (e.g., P=(I-bkgd)/c)


%PTC last updated 4/22/2020

nQ=length(Qx);

%integral(1/V exp(i Rtrans.q.r)) for a cylinder (radius=R, length=L)
%Qx,Qy,Qz q in cartesian
%order of rotations: theta, phi
cylinder=@(Qx,Qy,Qz,R,L,gamma,beta,alpha)4.*L.^(-1).*R.^(-1).*...
  besselj(1,R.*((Qy.*cos(alpha)+(-1).*Qx.*sin( ...
  alpha)).^2+(cos(beta).*(Qx.*cos(alpha)+Qy.*sin(alpha))+(-1).*Qz.* ...
  sin(beta)).^2).^(1/2)).*(Qz.*cos(beta)+(Qx.*cos(alpha)+Qy.*sin( ...
  alpha)).*sin(beta)).^(-1).*((Qy.*cos(alpha)+(-1).*Qx.*sin(alpha)) ...
  .^2+(cos(beta).*(Qx.*cos(alpha)+Qy.*sin(alpha))+(-1).*Qz.*sin( ...
  beta)).^2).^(-1/2).*sin((1/2).*L.*(Qz.*cos(beta)+(Qx.*cos(alpha)+ ...
  Qy.*sin(alpha)).*sin(beta)));


%Matrix initialization
P=zeros(nQ,1);
Rmat=P;
Lmat=P;
alphamat=P;
betamat=P;
gammamat=P;
Rmat(:)=R;
Lmat(:)=L;
cumulativeweight=0;
dangle=pi./nangle;

for Jp=1:nangle
    phi=(Jp-0.5)*dangle;
    for Jt=1:nangle
        theta=(Jt-0.5)*dangle;
        %setup matricies
        gammamat(:)=0;
        betamat(:)=theta;
        alphamat(:)=phi;
        cumulativeweight=cumulativeweight+2.*dangle.*dangle.*...
            (sin(theta)).*Y(Jt,Jp);
        %evaluate F at Qz=0 with weighting from simulation
        %flow-gradient
        P=P+2.*dangle.*dangle.*(sin(theta)).*Y(Jt,Jp).*...
            (cylinder(Qx,Qy,Qz,Rmat,Lmat,gammamat,betamat,alphamat)).^2;
    end
end

c=(10.^8).*pi.*R.^2.*L.*volfrac*(dsld)^2;
I=c.*P./cumulativeweight+bkgd;%calculate intensity in units of (cm^-1)

end

