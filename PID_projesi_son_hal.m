R=2.8;
La=1.4;
Ref=24;
E=Ref;
Kb=1.24;
K_i=1.24;
j=0.072;
B=0.02;
Ref_max=24;
Ref_min=24;

ACTF=tf(K_i,[La*j B*La+R*j B*R+K_i*Kb]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sureki zaman %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Basamak ess_s yi bulalım
pay=Ref;
payda=1+ACTF;
ess_s_basamak=pay/payda;
ess_deger_at_basamak=evalfr(ess_s_basamak,0);

disp(['sürekli zaman  ess_s basamak= ',num2str(ess_deger_at_basamak)]);

% rampa ess_s yi bulalım

s=tf('s');

pay=Ref/s;
payda=1+ACTF;
ess_s_rampa=pay/payda;
ess_deger_at_rampa=evalfr(ess_s_rampa,0);
disp(['sürekli zaman  ess_s_rampa= ',num2str(ess_deger_at_rampa)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T bul %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,d]=tfdata(ACTF);
n = cell2mat(n);  % 'n' hücresini sayısal vektöre dönüştürme n = {1, 2, 3};    % n bir hücre dizisi 
d = cell2mat(d);  % 'd' hücresini sayısal vektöre dönüştürme n = cell2mat(n);   % n şimdi bir sayısal vektördür: n = [1 2 3]
payda_kokleri = roots(d);
min_kok=max(abs(payda_kokleri));
kt=10;
T=1/(kt*min_kok);
Ts=T;
T=10e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ayrık zaman %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% z donusumu  ///Gp(z)
z_donusum_ACTF= c2d(ACTF,T,'zoh');
Gp=z_donusum_ACTF;
[nz,dz]=tfdata(z_donusum_ACTF);
nz= cell2mat(nz);  % 'n' hücresini sayısal vektöre dönüştürme n = {1, 2, 3};    % n bir hücre dizisi 
dz= cell2mat(dz);  % 'd' hücresini sayısal vektöre dönüştürme n = cell2mat(n);   % n şimdi bir sayısal vektördür: n = [1 2 3]


% basamak ess_z'i bulalım
z=tf('z',T);

pay=Ref*z;
payda=1+z_donusum_ACTF;
ess_z_basamak=pay/payda;
ACTF_Z=minreal(ess_z_basamak) ;

KCTF=simplify(ACTF_Z/(1+ACTF_Z));

ess_z_deger_at_basamak=evalfr(ess_z_basamak,1);

disp(['ayrık zaman    ess_z basamak= ',num2str(ess_z_deger_at_basamak)]);

% Rampa ess_z'i bulalım
z=tf('z',T);

pay=((Ref*z*T)/(z-1));
payda=1+z_donusum_ACTF;
ess_z_rampa=pay/payda;
ess_z_deger_at_rampa=evalfr(ess_z_rampa,1);

disp(['ayrık zaman    ess_z ayrık= ',num2str(ess_z_deger_at_rampa)]);

%%%%%%%%%%%%%%%%%%%%%%%%% ister %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yuzde_asim=4.32;   %yüz ile çarpılmış
a=pi^2+log(yuzde_asim/100)^2;
ksi=-log(yuzde_asim/100)/sqrt(a);  %ksi=0.707;
% b=-ksi*pi/sqrt(1-ksi^2);
% asim=100*exp(b);  %asim=4.3
ts=4/min_kok; % yüzde 2 kiriterine gore 4 tau zamanda yrleşsin
wn=4/(ksi*ts);


karakteristik_denk=[1 2*ksi*wn wn^2];
baskin_kok= roots(karakteristik_denk);
ayrik_zaman_baskin_kokleri= exp(baskin_kok * T);
birinci_kok=ayrik_zaman_baskin_kokleri(1);

Kh=1;   % ess yi küçültürsem daha hızlı sisteme oturacaktır büyültürsem ise daha yavaş
es_s=(Kh*2*ksi*Ref)/wn;  
es_z=es_s;

Kv=Ref/es_z;
limit=evalfr(z_donusum_ACTF,1);
Ki=Kv*T/limit;

z1=abs(birinci_kok);        %genlik
beta=angle(birinci_kok);   %açı

Gp_z=evalfr(Gp,birinci_kok);
Gp_z1=abs(Gp_z);  %genlik
fi=angle(Gp_z);   %acı


Kd=(z1/sin(beta))*(((Ki*sin(beta))/(z1-2*cos(beta)+(1/z1)))+(sin(fi)/Gp_z1));
Kp=(-cos(fi)/Gp_z1)-(2*Ki*z1*((z1-cos(beta))/(z1^2-2*z1*cos(beta)+1)))+((-z1*sin(fi)+cos(beta)*sin(fi))/(Gp_z1*sin(beta)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%% karalılık incelemesi%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z=tf('z',T);
KCTF_PID=((Kp+(Ki*z/(z-1))+Kd*(z-1)/z)*Gp)/(1+((Kp+(Ki*z/(z-1))+Kd*(z-1)/z)*Gp));
[pay, payda]=tfdata(KCTF_PID,'v');
pay_kokleri=roots(pay);
payda_kokleri=roots(payda);
ayrik_zaman_baskin_kokleri_V = [ayrik_zaman_baskin_kokleri(1), ayrik_zaman_baskin_kokleri(2)];  % Bu kısm wn^2s^2 +ksi Wns+

%%%%%%%%%%%%%%% Birim çemberin parametrelerini çizim için ayarlama %%%%%%%%%
theta = linspace(0, 2*pi, 100);
x = cos(theta);
y = sin(theta);

% Birim çemberin çizimi
figure;
plot(x, y, 'k--', 'LineWidth', 1); % Birim çember
hold on;

% Pay köklerini (real ve imaginary kısmını) çiz
plot(real(pay_kokleri), imag(pay_kokleri), 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Pay kökleri kırmızı renkte
disp('Payın Kökleri:');
disp(pay_kokleri);

% Payda köklerini (real ve imaginary kısmını) çiz
plot(real(payda_kokleri), imag(payda_kokleri), 'bx', 'MarkerSize', 8, 'LineWidth', 2); % Payda kökleri mavi renkte
disp('Paydanın Kökleri:');
disp(payda_kokleri);

% Ayrık zaman baskın köklerini (real ve imaginary kısmını) çiz
plot(real(ayrik_zaman_baskin_kokleri_V), imag(ayrik_zaman_baskin_kokleri_V), 'g*', 'MarkerSize', 8, 'LineWidth', 2); % Ayrık zaman baskın kökleri yeşil renkte
disp('Ayrık Zaman Baskın Kökleri_V:');
disp(ayrik_zaman_baskin_kokleri_V);

% Grafiği güzelleştirme
axis equal;
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
grid on;
title('Birim Çember Üzerinde Kökler');
xlabel('Reel');
ylabel('Sanal');
legend('Birim Çember', 'Pay Kökleri', 'Payda Kökleri', 'Ayrık Zaman Baskın Kökleri_V', 'Location', 'Best');
hold off;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Kp');
disp(Kp);
disp('Ki');
disp(Ki);
disp('Kd');
disp(Kd);