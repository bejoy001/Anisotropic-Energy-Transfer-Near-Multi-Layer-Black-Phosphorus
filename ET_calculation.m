clc;
clear all;
close all;
tic
e = 1.6*1e-19;
h = 6.634*1e-34;
hcut = h/(2*pi);
kB = 1.380649*1e-23;
T = 300;
kBT = kB*T;

% index=43;  %bulk
index= 43; % nl=60
Omega=logspace(13,17,100);
omega=Omega(index);

rD = [0, 0, 5*1e-9];
x0 = rD(1); y0 = rD(2); z0=rD(3);
x_plane=logspace(-9,-6,30);
y_plane=logspace(-9,-6,30);
% y_plane = 1e-9;
z_plane=5e-9;


eA = [0 0 1];
eD = [0 0 1];
eAx = eA(1);eAy = eA(2);eAz = eA(3);
eDx = eD(1);eDy = eD(2);eDz = eD(3);



m0 = 9.1*1e-31;
mu0 = 4*pi*1e-7;
epsilon0 = 8.854*1e-12;
c=1/sqrt(epsilon0*mu0);
eta0 = sqrt(mu0/epsilon0);
f = omega/(2*pi);
lambda = c./f;
k = omega/c;
k0 = (2*pi)./lambda;



K_max=1000;
N =10*1e3;
kmin = -K_max;
kmax = K_max;





load cond_xx.mat;
arysigxx=arysigxx(100,:);
load cond_yy.mat;
arysigyy=arysigyy(100,:);


% 
sigma_0=e^2/(4*hcut);


sigmaX=sigma_0*arysigxx(index);
sigmaY=sigma_0*arysigyy(index);

kx = linspace(kmin*k0, kmax*k0, N);
ky = linspace(kmin*k0, kmax*k0, N);

ET= zeros(length(y_plane),length(x_plane));
Gamma0=zeros(length(y_plane),length(x_plane));
Gamma=zeros(length(y_plane),length(x_plane));
ET_unnorm=zeros(length(y_plane),length(x_plane));
Vaccum=zeros(length(y_plane),length(x_plane));
PF=zeros(length(y_plane),length(x_plane));

for kk=1:length(y_plane)
    y = y_plane(kk);

    for j = 1:length(x_plane)

        fprintf('y = %d,x = %d Complete\n\n',kk,j);
        x = x_plane(j);
        z=z_plane;
        rA=[x y z];
        x=rA(1);y=rA(2); z=rA(3);

        Rmag=sqrt((x-x0).^2+(y-y0).^2+(z-z0).^2);
        Rvec = rD-rA;
        internal_int = zeros(1,N);

        for i = 1:length(kx)
            krho = sqrt(kx(i).^2+ky.^2);
            kz = sqrt(k0.^2-krho.^2);

            sigmaxx = (kx(i).^2*sigmaX+ky.^2*sigmaY)./krho.^2;
            sigmaxy = (kx(i).*ky.*(sigmaY-sigmaX))./krho.^2; %
            sigmayx = sigmaxy;
            sigmayy= (ky.^2*sigmaX+kx(i).^2*sigmaY)./krho.^2;
            Zs = kz./k0;
            Zp = k0./kz;
            cp = kz./k0;
            den_r = (2*Zs+eta0.*sigmayy).*(2*Zp+eta0.*sigmaxx)-eta0^2.*sigmaxy.*sigmayx;

            Rss = (-eta0.*sigmayy.*(2*Zp+eta0.*sigmaxx)+eta0.^2.*sigmaxy.*sigmayx)./den_r;
            Rsp = -(2.*cp.*Zp.*eta0.*sigmaxy)./den_r;
            Rps = (2.*Zs.*eta0.*sigmayx)./(cp.*den_r);
            Rpp = (eta0.*sigmaxx.*(2.*Zs+eta0.*sigmayy)+eta0.^2.*sigmaxy.*sigmayx)./den_r;

            %     M11 = (Rss.*(ky.^2))./(kz.*krho.^2)+((-kx(i).*ky).*Rsp)./(k0.*krho.^2)+((kx(i).*ky).*Rps)./(k0.*krho.^2)-(Rpp.*(kx(i)^2.*kz))./(k0.^2.*krho.^2);
            %     M12 = (Rss.*(-kx(i).*ky))./(kz.*krho.^2)+((-ky.^2).*Rsp)./(k0.*krho.^2)+((-kx(i).^2).*Rps)./(k0.*krho.^2)-(Rpp.*(kx(i).*ky.*kz))./(k0^2.*krho.^2);
            %     M13 = (((-ky.*krho.^2)./kz).*Rsp)./(k0.*krho.^2)-(Rpp.*((kx(i).*krho.^2.*kz)./kz))./(k0^2.*krho.^2);
            %     M21 = (Rss.*(-kx(i).*ky))./(kz.*krho.^2)+((kx(i).^2).*Rsp)./(k0.*krho.^2)+((ky.^2).*Rps)./(k0.*krho.^2)-(Rpp.*(kx(i).*ky.*kz))./(k0^2.*krho.^2);
            %     M22 = (Rss.*(kx(i).^2))./(kz.*krho.^2)+((kx(i).*ky).*Rsp)./(k0.*krho.^2)-((kx(i).*ky).*Rps)./(k0.*krho.^2)-(Rpp.*(ky.^2.*kz))./(k0^2.*krho.^2);
            %     M23 = (((kx(i).*krho.^2)./kz).*Rsp)./(k0.*krho.^2)-(Rpp.*((ky.*krho.^2)))./(k0^2.*krho.^2);
            %     M31 = -(((ky.*krho.^2)./kz).*Rps)./(k0.*krho.^2)+(Rpp.*(kx(i).*krho.^2))./(k0^2.*krho.^2);
            %     M32 = (((kx(i).*krho.^2)./kz).*Rps)./(k0.*krho.^2)+(Rpp.*(ky.*krho.^2))./(k0^2.*krho.^2);
            M33 = (Rpp.*(krho.^2))./(k0^2.*kz);
            fun =M33;
            integrand = fun.*(exp(1i*kx(i)*(x-x0)+1i*ky*(y-y0))).*exp(1i*kz*(z+z0));
            internal_int(i) = trapz(ky, integrand);

            integrand_PF = fun.*exp(1i*kz*(z0+z0));
            internal_int_PF(i) = trapz(ky, integrand_PF);

        end
        subplot(2,1,1);semilogy(kx/k0,abs(internal_int));
        subplot(2,1,2);semilogy(kx/k0,abs(internal_int_PF));
        Gamma(kk,j) = 1i*trapz(kx, internal_int)/(8*pi^2);
        %         Gamma0(kk,j) = exp(1i*k0.*Rmag)./(4*pi*Rmag).*(1+(1i*k0.*Rmag-1)./(k0.^2.*Rmag.^2));
        eA=[0 0 1];eD=eA;
        term1 = (1+(1i*k0.*Rmag-1)./(k0.^2.*Rmag.^2))*eye(3);
        term2 = (3-3*1i*k0.*Rmag-k0.^2.*Rmag.^2)./(k0.^2.*Rmag.^2)*(Rvec'*Rvec)/(Rmag.^2);
        Gvac = (exp(1i*k0.*Rmag)./(4*pi*Rmag))*(term1+term2);
        Green_Vac=eA*(Gvac*eD');
        Gamma0(kk,j) = Green_Vac;
        ET_unnorm(kk,j) = abs(Gamma(kk,j))^2;
        Vaccum(kk,j)= abs(Gamma0(kk,j))^2;
        ET(kk,j) = abs(Gamma(kk,j)+Gamma0(kk,j))^2/abs(Gamma0(kk,j))^2;
        Green_PF=1i*trapz(kx, internal_int_PF)/(8*pi^2);
        PF(kk,j)=1+(6*pi*c./omega).*imag(Green_PF);
    end
end


figure
subplot(4,1,1);loglog(x_plane/1e-9,ET_unnorm,'Linewidth',1.5);
hold on;loglog(x_plane/1e-9,Vaccum,'Linewidth',1.5);legend('scattered','vaccum');% title('rA= [ x y 0] ;   rD=[0 0 5] nm;; y=10 nm');
grid on

subplot(4,1,2);loglog(x_plane/1e-9,ET_unnorm,'Linewidth',2);
% legend(' Scattered, z= 1 nm','z=10 nm', 'z=100 nm')
legend(' Scattered, y= 1 nm','y=10 nm', 'y=100 nm')

subplot(4,1,3)
semilogx(x_plane/1e-9,ET,'Linewidth',2);legend('Scattered/Vacum');xlabel('x (nm)');grid on
% legend('z= 1 nm','z=10 nm','z=100 nm','z=1000 nm')
legend('y= 1 nm' ,'y=10 nm', 'y=100 nm')

subplot(4,1,4);semilogy(kx/k0,abs(internal_int))

toc