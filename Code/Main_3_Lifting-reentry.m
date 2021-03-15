%% PROGRAMA 3 - CAS DE REENTRADA SUSTENTADORA

% INSTRUCCIONS:
%     1. Escollir les gràfiques a mostrar per pantalla introduint "false" o "true" al vector "grafiques_a_mostrar"
%     2. Indicar amb un "false" o "true" a la variable "printear" si es vol que les gràfiques mostrades es guardin automàticament com a .png
%     3. Es poden modificar els valors del Flight Path Angle (gamma) a tractar. Si es volen afegir més valors, afegir més components al vector "gamma_0_vector_deg"

clear
clc
close all
format long

grafiques_a_mostrar=[true;  % (1): Velocitat
                    true;   % (2): Temps d'entrada
                    true;   % (3): Abast
                    true;   % (4): Desacceleració
                    true;   % (5): Pressió dinàmica
                    true;   % (6): Mach
                    true;   % (7): Reynolds
                    true;   % (8): Pressió estancament
                    true;   % (9): Entalpia estancament
                    true;   % (10): Flux de calor al punt d'estancament
                    true;   % (11): Energia dinàmica
                    false;   % (12): Angle de descens
                    false];   % (13): Gràfiques de les propietats de l'atmosfera (T(h),densitat(h),etc)

printear=false;

% Vector de gammes avaluades
gamma_0_vector_deg=[0.1 1 2.5]; % [graus]

%% Càlculs previs

% Definició de constants físiques
R_asterisc=8.31432; % [N m/(mol K)]
M_0=28.9644*1e-3; % [kg/mol] 
N_A=6.022169e26*1e-3; % [1/mol]
g_0_geo=9.80665; % [m^2/(s^2 m']
g_0=9.80665; % [m/s^2]
Gamma=g_0_geo/g_0; % [m'/m]
r_0=6.356766e6; % [m]
H_b=[0 11 20 32 47 51 71 84.8520]*1e3; % [m']
L_Mb=[-6.5 0.0 1.0 2.8 0.0 -2.8 -2.0]*1e-3; % [K/m']
beta_viscosa=1.458*10^(-6); % [kg /(s m K^1/2]
S=110.4; %[K] Sutherland's constant
gamma_aire=1.4;
eficiencia=1; % (=C_L/C_D)
L=107.5/unitsratio('feet', 'meter'); % [m]

% Càlcul de pressions base (P_b) i temperatures base a escala molecular (T_b) de cada capa de l'atmosfera
% en base a U.S. Standard Atmosphere 1976
T_Mb=[288.15 zeros(1,length(L_Mb)-1)]; % [K]
P_b=[101325.0 zeros(1,length(L_Mb)-1)]; % [Pa]
for i=1:(length(L_Mb)-1)
    T_Mb(i+1)=T_Mb(i)+L_Mb(i)*(H_b(i+1)-H_b(i));
    if L_Mb(i)==0
        P_b(i+1)=P_b(i)*exp((-g_0_geo*M_0*(H_b(i+1)-H_b(i)))/(R_asterisc*T_Mb(i)));
    else
        P_b(i+1)=P_b(i)*(T_Mb(i)/(T_Mb(i)+L_Mb(i)*(H_b(i+1)-H_b(i))))^(g_0_geo*M_0/(R_asterisc*L_Mb(i)));
    end
end

% Condicions inicials
z_0=250000/unitsratio('feet', 'meter');
V_0=23000/unitsratio('feet', 'meter');
range_0=0;
t_0=0;

% Definició de l'interval de la variable independent (variació del temps)
delta_t=0.25;

% Definició del vector de la variable independent
t=t_0:delta_t:1500;

% Es determina el nombre d'elements de les matrius de resultats
N=length(t);

% Inicialització de matrius d'emmagatzemament de dades
gamma=zeros(1,N);
V=[V_0 zeros(1,N-1)];
h=[z_0 zeros(1,N-1)];
range=[range_0 zeros(1,N-1)];

%Definició vector d'altura, per obtenció d'informació.
N_Z=100000;
Z=linspace(84851,0,N_Z);

% Càlcul de l'altitud geopotencial H segons Z
H=Z*Gamma*r_0./(r_0+Z);

% Inicialització de matrius d'emmagatzemament de dades
P=zeros(1,N_Z);
T_M=zeros(1,N_Z);
P_t=zeros(1,N);

% Càlcul de temperatures i pressions segons l'altura geopotencial H
% en base a U.S. Standard Atmosphere 1976
for i=1:N_Z
for b=1:(length(L_Mb))
    if H(i)>=H_b(b) && H(i)<H_b(b+1)
       T_M(i)=T_Mb(b)+L_Mb(b)*(H(i)-H_b(b));
       if L_Mb(b)==0
           P(i)=P_b(b)*exp((-g_0_geo*M_0*(H(i)-H_b(b)))/(R_asterisc*T_Mb(b)));
       else
           P(i)=P_b(b)*(T_Mb(b)/(T_Mb(b)+L_Mb(b)*(H(i)-H_b(b))))^(g_0_geo*M_0/(R_asterisc*L_Mb(b)));
       end
       break
    end
end
end

% Càlcul de la gravetat segons Z
g=g_0*(r_0./(r_0+Z)).^2;

% Càlcul d'altres propietats termofísiques
T=T_M; % Fins a 80 km
rho=(P*M_0)./(T*R_asterisc); 
mu=(beta_viscosa*T.^(1.5))./(S+T); 

% Càlcul de la velocitat del so segons T_M
a=sqrt(gamma_aire*R_asterisc*T_M/M_0);

% Inicialització dels vectors amb els valors particulars de les propietats
% termofísiques per cada posicio en la trajectoria

M_particular=zeros(1,N);
Temperatura_particular=zeros(1,N);
rho_particular=zeros(1,N);
mu_particular=zeros(1,N);
P_particular=zeros(1,N);
a_particular=zeros(1,N);

%% SOLVER PER RUNGE-KUTTA DE QUART ORDRE

% Conversió d'unitats de la beta
beta=2e5/(2690*0.84)*unitsratio('feet', 'meter')^2*4.4482216152605;

%Conversió d'unitats de la gamma
gamma_0_vector=deg2rad(gamma_0_vector_deg);

for j=1:length(gamma_0_vector)
    gamma(1)=gamma_0_vector(j);

% Declaració de les equacions diferencials
F_Vt=@(t_v,V_v,gamma_v,h_v,range_v,g_v,a_v,P_v) g_v*(-(0.5*( V_v/a_v )^2*gamma_aire*P_v )/beta+sin(gamma_v));
F_gammat=@(t_v,V_v,gamma_v,h_v,range_v,g_v,a_v,P_v) (-(0.5*( V_v/a_v )^2*gamma_aire*P_v )*g_v/beta*eficiencia+ cos(gamma_v)*(g_v-(V_v)^2/(r_0+h_v)))/V_v;
F_ht=@(t_v,V_v,gamma_v,h_v,range_v) -V_v*sin(gamma_v);
F_ranget=@(t_v,V_v,gamma_v,h_v,range_v) r_0*V_v*cos(gamma_v)/(r_0+h_v);

for i=1:N-1 % Bucle de resolució per Runge-Kutta
    
    % Cerca de les propietats termofísiques associades a la posició actual
    pos=find((abs(Z-h(i)))==min(abs(Z-h(i))));
    
    % Emmagatzematge de les propietats termofísiques puntuals
    a_particular(i)=a(pos);
    Temperatura_particular(i)=T(pos);
    rho_particular(i)=rho(pos);
    mu_particular(i)=mu(pos);
    P_particular(i)=P(pos);
        
    % Coeficients 1
    k_1Vt =     F_Vt(t(1,i),    V(1,i),gamma(1,i),h(1,i),range(1,i),g(pos),a(pos),P(pos));
    k_1gammat = F_gammat(t(1,i),V(1,i),gamma(1,i),h(1,i),range(1,i),g(pos),a(pos),P(pos));
    k_1ht =     F_ht(t(1,i),    V(1,i),gamma(1,i),h(1,i),range(1,i));
    k_1ranget = F_ranget(t(1,i),V(1,i),gamma(1,i),h(1,i),range(1,i));
    
    %Coeficients 2
    k_2Vt=      F_Vt(t(1,i)+delta_t/2,      V(1,i)+delta_t/2*k_1Vt, gamma(1,i)+delta_t/2*k_1gammat, h(1,i)+delta_t/2*k_1ht, range(1,i)+delta_t/2*k_1ranget,g(pos),a(pos),P(pos));
    k_2gammat=  F_gammat(t(1,i)+delta_t/2,  V(1,i)+delta_t/2*k_1Vt, gamma(1,i)+delta_t/2*k_1gammat, h(1,i)+delta_t/2*k_1ht, range(1,i)+delta_t/2*k_1ranget,g(pos),a(pos),P(pos));
    k_2ht=      F_ht(t(1,i)+delta_t/2,      V(1,i)+delta_t/2*k_1Vt, gamma(1,i)+delta_t/2*k_1gammat, h(1,i)+delta_t/2*k_1ht, range(1,i)+delta_t/2*k_1ranget);
    k_2ranget=  F_ranget(t(1,i)+delta_t/2,  V(1,i)+delta_t/2*k_1Vt, gamma(1,i)+delta_t/2*k_1gammat, h(1,i)+delta_t/2*k_1ht, range(1,i)+delta_t/2*k_1ranget);
    
    % Coeficients 3
    k_3Vt=      F_Vt(t(1,i)+delta_t/2,      V(1,i)+delta_t/2*k_2Vt, gamma(1,i)+delta_t/2*k_2gammat, h(1,i)+delta_t/2*k_2ht, range(1,i)+delta_t/2*k_2ranget,g(pos),a(pos),P(pos));
    k_3gammat=  F_gammat(t(1,i)+delta_t/2,  V(1,i)+delta_t/2*k_2Vt, gamma(1,i)+delta_t/2*k_2gammat, h(1,i)+delta_t/2*k_2ht, range(1,i)+delta_t/2*k_2ranget,g(pos),a(pos),P(pos));
    k_3ht=      F_ht(t(1,i)+delta_t/2,      V(1,i)+delta_t/2*k_2Vt, gamma(1,i)+delta_t/2*k_2gammat, h(1,i)+delta_t/2*k_2ht, range(1,i)+delta_t/2*k_2ranget);
    k_3ranget=  F_ranget(t(1,i)+delta_t/2,  V(1,i)+delta_t/2*k_2Vt, gamma(1,i)+delta_t/2*k_2gammat, h(1,i)+delta_t/2*k_2ht, range(1,i)+delta_t/2*k_2ranget);
    
    % Coeficients 4
    k_4Vt=      F_Vt(t(1,i)+delta_t,        V(1,i)+delta_t*k_3Vt,   gamma(1,i)+delta_t*k_3gammat,   h(1,i)+delta_t*k_3ht,   range(1,i)+delta_t*k_3ranget,g(pos),a(pos),P(pos));
    k_4gammat=  F_gammat(t(1,i)+delta_t,    V(1,i)+delta_t*k_3Vt,   gamma(1,i)+delta_t*k_3gammat,   h(1,i)+delta_t*k_3ht,   range(1,i)+delta_t*k_3ranget,g(pos),a(pos),P(pos));
    k_4ht=      F_ht(t(1,i)+delta_t,        V(1,i)+delta_t*k_3Vt,   gamma(1,i)+delta_t*k_3gammat,   h(1,i)+delta_t*k_3ht,   range(1,i)+delta_t*k_3ranget);
    k_4ranget=  F_ranget(t(1,i)+delta_t,    V(1,i)+delta_t*k_3Vt,   gamma(1,i)+delta_t*k_3gammat,   h(1,i)+delta_t*k_3ht,   range(1,i)+delta_t*k_3ranget);
    
    % Cálcul dels valors posteriors
    V(1,i+1)=       V(1,i)+     delta_t/6*(k_1Vt+2*k_2Vt+2*k_3Vt+k_4Vt);
    gamma(1,i+1)=   gamma(1,i)+ delta_t/6*(k_1gammat+2*k_2gammat+2*k_3gammat+k_4gammat);
    h(1,i+1)=       h(1,i)+     delta_t/6*(k_1ht+2*k_2ht+2*k_3ht+k_4ht);
    range(1,i+1)=   range(1,i)+ delta_t/6*(k_1ranget+2*k_2ranget+2*k_3ranget+k_4ranget);
end

% Càlcul del nombre de Mach
M=V./a_particular;

% Càlcul de la pressió dinàmica
Q=0.5*M.^2*gamma_aire.*P_particular;

% Càlcul del nombre de Reynolds
Re=V.*rho_particular*L./mu_particular;

% Càlcul de l'acceleració
V_punt=-Q./beta+sin(gamma);

% Càlcul de l'energia dinàmica
energia_dinamica=Q.*V;

for i=1:N % Bucle per a les pressions d'estancament

if (M(i)>1)%sqrt((gamma_aire-1)/(2*gamma_aire)))  % Règim supersònic (ona de xoc)  
p_shock=P_particular(i)*(1+((2*gamma_aire)/(gamma_aire+1))*(M(i)^2-1));
P_t(i)=p_shock*(1+((gamma_aire-1)/(2))*((1+((gamma_aire-1)/2)*M(i)^2)/(gamma_aire*M(i)^2-((gamma_aire-1)/2)))).^(gamma_aire/(gamma_aire-1));
else % Règim subsònic
    P_t(i)=P_particular(i)*(1+(gamma_aire-1)*0.5*M(i)^2).^(gamma_aire/(gamma_aire-1));
end % Fi del condicional

end % Fi del bucle de les pressions d'estancament

% Fórmula del C_p del Formulae
C_p=1034.9-2.849*10^(-1)*Temperatura_particular+7.817*10^(-4)*Temperatura_particular.^2-4.971*10^(-7)*Temperatura_particular.^3+1.007*10^(-10)*Temperatura_particular.^4;

% Càlcul de l'entalpia d'estancament
entalpia_estancament=C_p.*Temperatura_particular+0.5* V.^2;

% Flux de calor per la correlació de Sutton and Graves
q_SG=1.83*10^(-4)*sqrt(rho_particular./(1/unitsratio('feet', 'meter'))).*V.^3;


%% IMPRESSIÓ DE GRÀFIQUES

if grafiques_a_mostrar(1)==true
fig1=figure(1); % Velocitat
set(fig1,'Renderer', 'painters', 'Position', [15 565 375 170]);
plot(V,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Velocitat $V\;\mathrm{(m/s)}$','Interpreter','latex')
legend('Location','southeast','Interpreter','latex')

if(j==1)
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
end

hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f1','velocitat_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(2)==true
fig2=figure(2); % Temps d'entrada
set(fig2,'Renderer', 'painters', 'Position', [400 565 375 170]);
plot(t,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Temps d''entrada $t\;\mathrm{(s)}$','Interpreter','latex')
legend('Location','southwest','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f2','temps_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(3)==true
fig3=figure(3);  % Abast
set(fig3,'Renderer', 'painters', 'Position', [785 565 375 170]);
plot(range/1000,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Abast $r\;\mathrm{(km)}$','Interpreter','latex')
legend('Location','southwest','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f3','abast_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(4)==true
fig4=figure(4); % Desacceleració
set(fig4,'Renderer', 'painters', 'Position', [1170 565 375 170]);
plot(-V_punt,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Desacceleraci\''o $-\dot{V}\;\mathrm{(g)}$','Interpreter','latex')
legend('Location','southeast','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f4','desacceleracio_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(5)==true
fig5=figure(5); % Pressió dinàmica
set(fig5,'Renderer', 'painters', 'Position', [15 305 375 170]);
plot(Q,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Pressi\''o din\`amica $Q\;\mathrm{(N/{m}^2)}$','Interpreter','latex')
legend('Location','southwest','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f5','pressio_dinamica_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(6)==true
fig6=figure(6); % Mach
set(fig6,'Renderer', 'painters', 'Position', [400 305 375 170]);
plot(M,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('N\''umero de Mach $M$','Interpreter','latex')
legend('Location','southeast','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f6','mach_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(7)==true
fig7=figure(7); % Reynolds
set(fig7,'Renderer', 'painters', 'Position', [785 305 375 170]);
plot(Re,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Reynolds $Re$','Interpreter','latex')
legend('Location','southwest','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f7','reynolds_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(8)==true
fig8=figure(8); % Pressió estancament
set(fig8,'Renderer', 'painters', 'Position', [1170 305 375 170]);
plot(P_t,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Pressi\''o d''estancament $P_{t}\;\mathrm{(N/{m}^2)}$','Interpreter','latex')
legend('Location','northeast','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f8','pressio_estancament_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(9)==true
fig9=figure(9); % Entalpia estancament
set(fig9,'Renderer', 'painters', 'Position', [15 45 375 170]);
plot(entalpia_estancament,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Entalpia d''estancament $h_{t}\;\mathrm{(J/kg)}$','Interpreter','latex')
legend('Location','southeast','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f9','entalpia_estancament_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(10)==true
fig10=figure(10); % Flux de calor al punt d'estancament (segons Sutton Graves)
set(fig10,'Renderer', 'painters', 'Position', [400 45 375 170]);
plot(-q_SG,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Flux calor pel punt d''estancament $\dot{q}_{t}\;\mathrm{(J/m^2s)}$','Interpreter','latex')
legend('Location','southwest','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f10','flux_calor_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(11)==true
fig11=figure(11); % Energia dinàmica
set(fig11,'Renderer', 'painters', 'Position', [785 45 375 170]);
plot(energia_dinamica,h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$'])
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Energia din\`amica $e_{d}\;\mathrm{(J/m^2s)}$','Interpreter','latex')
legend('Location','southeast','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f11','energia_dinamica_lifting.png','-dpng','-r500');    
end
end

if grafiques_a_mostrar(12)==true
fig12=figure(12); % Angle de lliscament
set(fig12,'Renderer', 'painters', 'Position', [1170 45 375 170]);
plot(rad2deg(gamma),h/1000,'DisplayName',['$\gamma=$ ', num2str(gamma_0_vector_deg(j)),'$^\circ$']);
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Flight Path Angle $\gamma\;\mathrm{(^{\circ})}$','Interpreter','latex')
xticks([0,15,30,45,60,75,90]);
legend('Location','northeast','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on
if(j==length(gamma_0_vector) && printear==true)
print('-f12','gamma_lifting.png','-dpng','-r500');    
end
end

end % Fi del for del solver

% Impressió de gràfiques de l'evolució de propietats termofísiques amb
if grafiques_a_mostrar(13)==true
fig13=figure(13); % temperatura en funció de l'alçada
set(fig13,'Renderer', 'painters', 'Position', [500 150 550 550]);
plot(T,Z/1000,'r','linewidth',1.25)
title('Evoluci\''o de la temperatura amb l''altura','interpreter','latex');
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Temperatura $\mathrm{(K)}$','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on

fig14=figure(14); % pressió en funció de l'alçada
set(fig14,'Renderer', 'painters', 'Position', [500 150 550 550]);
plot(P,Z/1000,'r','linewidth',1.25)
title('Evoluci\''o de la pressi\''o amb l''altura','interpreter','latex');
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Pressi\''o $\mathrm{(Pa)}$','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on

fig15=figure(15); % pressió en funció de l'alçada
set(fig15,'Renderer', 'painters', 'Position', [500 150 550 550]);
plot(a,Z/1000,'r','linewidth',1.25)
title('Evoluci\''o de la velocitat del so amb l''altura','interpreter','latex');
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Velocitat del so $\mathrm{(m/s)}$','Interpreter','latex')
%legend('Location','southeast','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on

fig16=figure(16); % pressió en funció de l'alçada
set(fig16,'Renderer', 'painters', 'Position', [500 150 550 550]);
plot(rho,Z/1000,'r','linewidth',1.25)
title('Evoluci\''o de la densitat amb l''altura','interpreter','latex');
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Densitat $\mathrm{(kg/m^3)}$','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on

fig17=figure(17); % pressió en funció de l'alçada
set(fig17,'Renderer', 'painters', 'Position', [500 150 550 550]);
plot(mu,Z/1000,'r','linewidth',1.25)
title('Evoluci\''o de la viscositat din\`amica amb l''altura','interpreter','latex');
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Viscositat din\`amica $\mathrm{(kg/m\,s)}$','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on

C_p_gen=1034.9-2.849*10^(-1)*T+7.817*10^(-4)*T.^2-4.971*10^(-7)*T.^3+1.007*10^(-10)*T.^4;
fig18=figure(18); % calor a pressió constant en funció de l'alçada
set(fig18,'Renderer', 'painters', 'Position', [500 150 550 550]);
plot(C_p_gen,Z/1000,'r','linewidth',1.25)
title('Evoluci\''o de la calor espec\''ifica a pressi\''o constant amb l''altura','interpreter','latex');
ylabel('Altitud $h\;\mathrm{(km)}$','Interpreter','latex')
xlabel('Calor espec\''ifica a pressi\''o constant $\mathrm{(J/kg\,K)}$','Interpreter','latex')
ylim([0 85])
grid on
grid minor
ax = gca;
ax.GridColor = [0, 0, 0];
ax.GridAlpha=0.3;
ax.MinorGridColor = [0, 0, 0];
ax.MinorGridAlpha=0.3;
hold on

if (printear==true)   
    print('-f13','temperatura.png','-dpng','-r600');    
    print('-f14','pressio.png','-dpng','-r600');    
    print('-f15','velocitat_del_so.png','-dpng','-r600'); 
    print('-f16','densitat.png','-dpng','-r600');    
    print('-f17','viscositat.png','-dpng','-r600');    
    print('-f18','Cp.png','-dpng','-r600');
end
end % Fi del condicional de les gràfiques de l'evolució de propietats termofísiques

