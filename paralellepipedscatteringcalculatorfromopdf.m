function [I,c] = paralellepipedscatteringcalculatorfromopdf(Qx,Qy,Qz,Y,nangle,ngamma,a,b,c,volfrac,dsld,bkgd)
%%Function solves for scattering from dilute oriented paralellepipeds
%Run this after generating orientation distribution (Y)
%INPUTS
%Qx,Qy,Qz: arrays of the Qx,Qy,Qz coordinates where the scattering 
%    should be calculated in Angstroms^-1
%Y: The OPDF as a matrix generated from the FP simulations
%nangle: number of angles used as input to FP simulations
%ngamma: number of angles of gamma  to sample, the OPDF in this direction  
%   is always considered to be uniform in this implementation
%a: shortest length of paralellepiped in inverse units of Q (Angstroms)
%b: middle length of paralellepiped in inverse units of Q (Angstroms)
%c: longest length of paralellepiped in inverse units of Q (Angstroms)
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

parallelepiped= @(Qx,Qy,Qz,a,b,c,alpha,beta,gamma)8.*a.^(-1).*...
  b.^(-1).*c.^(-1).*(Qz.*cos(beta)+(Qx.*cos(alpha)+Qy.* ...
  sin(alpha)).*sin(beta)).^(-1).*(cos(gamma).*(cos(beta).*(Qx.*cos( ...
  alpha)+Qy.*sin(alpha))+(-1).*Qz.*sin(beta))+(Qy.*cos(alpha)+(-1).* ...
  Qx.*sin(alpha)).*sin(gamma)).^(-1).*(cos(gamma).*(Qy.*cos(alpha)+( ...
  -1).*Qx.*sin(alpha))+(-1).*(cos(beta).*(Qx.*cos(alpha)+Qy.*sin( ...
  alpha))+(-1).*Qz.*sin(beta)).*sin(gamma)).^(-1).*sin((1/2).*c.*( ...
  Qz.*cos(beta)+(Qx.*cos(alpha)+Qy.*sin(alpha)).*sin(beta))).*sin(( ...
  1/2).*a.*(cos(gamma).*(cos(beta).*(Qx.*cos(alpha)+Qy.*sin(alpha))+ ...
  (-1).*Qz.*sin(beta))+(Qy.*cos(alpha)+(-1).*Qx.*sin(alpha)).*sin( ...
  gamma))).*sin((1/2).*b.*(cos(gamma).*(Qy.*cos(alpha)+(-1).*Qx.* ...
  sin(alpha))+(-1).*(cos(beta).*(Qx.*cos(alpha)+Qy.*sin(alpha))+(-1) ...
  .*Qz.*sin(beta)).*sin(gamma)));
% integral(1/V exp(i Rtrans.q.r)) for a parallelepiped

amat=zeros(nQ,1);
amat(:)=a;
bmat=zeros(nQ,1);
bmat(:)=b;
cmat=zeros(nQ,1);
cmat(:)=c;
P=zeros(nQ,1);
gammamat=P;
betamat=P;
alphamat=P;
cumulativeweight=0;
dangle=pi./nangle;
dgamma=pi./ngamma;

for Jg=1:ngamma
    gamma=(Jg-0.5)*dgamma;
    for Jp=1:nangle
        phi=(Jp-0.5)*dangle;
        for Jt=1:nangle
            theta=(Jt-0.5)*nangle;
            %setup matricies
            gammamat(:)=gamma;
            betamat(:)=theta;
            alphamat(:)=phi;
            cumulativeweight=cumulativeweight+...
                2.*dangle.*dangle.*(sin(theta)).*Y(Jt,Jp)./ngamma;
            %evaluate F at Qz=0 with weighting from simulation
            %flow-gradient
            P=P+2.*dangle.*dangle.*(sin(theta)).*Y(Jt,Jp).*...
                (parallelepiped(Qx,Qy,Qz,amat,bmat,cmat,...
                alphamat,betamat,gammamat)).^2./ngamma;

        end
    end
end

c=(10.^8).*a.*b.*c.*volfrac*(dsld)^2;
I=c.*P./cumulativeweight+bkgd;%calculate intensity in units of (cm^-1)

end

