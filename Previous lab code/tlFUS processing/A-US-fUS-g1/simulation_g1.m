tau = (1:100)*0.2*1e-3;%s
vx = 5*1e-3;% m/s
vz = 10*1e-3;% m/s
k0 = 2*pi*18e6/1540; % Hz/(m/s)
p = 0.3;
thegmax = 50e-6; %m
thegmaz = 50e-6; %m
F = 1;

g1 = @(vx,vz) F*exp(-(vx*tau).^2/(4*thegmax^2)-(vz*tau).^2/(4*thegmaz^2)-(p*vz*k0*tau).^2).*exp(-i*2*k0*vz*tau);
m = 1; n = 1;clear g1s;
vx = [5]*1e-3;
vz = [5,10]*1e-3;
for m = 1: length(vx)
    for n = 1:length(vz)
g1s(m*n,:) = g1(vx(m),vz(n));
    end
end
absg1 = abs(g1s');
rg1 = real(g1s');
ig1 = imag(g1s');
figure; hold on;subplot(221);plot(sqrt(-log(absg1)));title('abs(g1)');
 hold on;subplot(222);plot(rg1);title('real(g1)');
hold on; subplot(223);plot(g1s');title('complex(g1)');
 hold on;subplot(224);plot(ig1);title('imag(g1)');
 hold off
 
 fig = figure;
 set(fig,'Position',[300,300,250,600])
 subplot(311);plot(tau*1e3,absg1(:,1),'-b','LineWidth',1);xlabel('\tau (ms)');ylabel('abs(g_{1})');
hold on;plot(tau*1e3,absg1(:,2),'-r','LineWidth',1);
 %  legend({'v_{x}=5mm/s, v_{z}=5mm/s','v_{x}=10mm/s, v_{z}=5mm/s'},'Location','bestoutside','Orientation','horizontal');
 legend({'v_{x}=5,v_{z}=5','v_{x}=5,v_{z}=10'});
 legend('boxoff');
 hold on;subplot(312);plot(tau*1e3,rg1(:,1),'-b','LineWidth',1);xlabel('\tau (ms)');ylabel('Re(g_{1})');
 hold on;plot(tau*1e3,rg1(:,2),'-r','LineWidth',1);
 hold on;subplot(313);plot(tau*1e3,ig1(:,1),'-b','LineWidth',1);xlabel('\tau (ms)');ylabel('Im(g_{1})');
 hold on;plot(tau*1e3,ig1(:,2),'-r','LineWidth',1);
 hold off
 
  fig = figure;
 set(fig,'Position',[300,300,250,250]);
 plot(g1s(1,:),'-k','LineWidth',1);xlabel('Re(g_{1})');ylabel('Im(g_{1})');
 grid on;