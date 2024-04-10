%% Inertia Script for PSET 1
clear; clc;

%% Center of Mass
% Bus
m_prism = 4.7578;
CoM_prism = [0 0 1.0019];
m_AV = 90.254;
CoM_AV = [0 0 0.204];
m_solar = 6*6.96;
CoM_solar =  [0 0 0.01];
m_beam = 67.23;
CoM_beam =  [0 0 0.9945];
m_bus = m_AV + m_prism + m_solar + m_beam;
CoM_bus = (m_AV * CoM_AV + m_beam * CoM_beam + m_solar * CoM_solar + m_prism * CoM_prism)/m_bus;

%Payload
m_I1 = 18;
CoM_I1 =  [-0.1375 0.1125 0.9945];
m_I2 = 25.4;
CoM_I2 = [0.25 0.175 0.9945];
m_I3 = 20.6;
CoM_I3 =  [0.205 -0.125 0.9945];
m_I4 = 22;
CoM_I4 = [-0.1975 -0.1525 0.9945];
m_payload = m_I1 + m_I2 + m_I3 + m_I4;
CoM_payload = (m_I1 *CoM_I1+m_I2 * CoM_I2 + m_I3 *CoM_I3 + m_I4 * CoM_I4)/m_payload;

%Total
CoM = (m_payload*CoM_payload + m_bus * CoM_bus)/(m_payload+m_bus);


%% Inertia

%prism
crosssec_prism = 3*sqrt(3)/2 *1^2;
h_prism = 1.2023;
I_y_prism = 5*sqrt(3)* 1^4/16; % = I_x_prism
rho_prism = m_prism/(crosssec_prism  * h_prism);
Ixx_prism = rho_prism * (h_prism * I_y_prism + crosssec_prism * h_prism^3/12);
I_zz_prism = 2*rho_prism * h_prism * I_y_prism;
I_prism = diag([Ixx_prism, Ixx_prism, I_zz_prism]);

%AV
h_AV = 0.4008;
rho_AV = m_AV/(crosssec_prism  * h_AV);
Ixx_AV = rho_AV * (h_AV * I_y_prism + crosssec_prism * h_AV^3/12);
I_zz_AV = 2*rho_AV * h_AV * I_y_prism;
I_AV = diag([Ixx_AV, Ixx_AV, I_zz_AV]);

%solar
I_solar = 1/12*m_solar/6 * diag([(0.5774^2 + 0.01^2), (1.1965^2 + 0.01^2), (0.5774^2 + 1.1965^2)]);

%IBeam
m_beam1 = 14.76;
m_beam2 = 37.3406;
m_beam3 = 14.76;

I_beam1 = 1/12*m_beam1 * diag([(0.535^2 + 1.1873^2), (0.03^2 + 1.1873^2), (0.535^2 + 0.03^2)]);
I_beam2 = 1/12*m_beam1 * diag([(0.05^2 + 1.1873^2), (0.85^2 + 1.1873^2), (0.05^2 + 0.85^2)]);
I_beam3 = I_beam1;

CoM_beam1 = [-0.44 0 0.9945];
CoM_beam2 = [0 0 0.9945];
CoM_beam3 = [0.44 0 0.9945];

%Payload
I_I1 = 1/12*m_I1 * diag([(0.175^2 + 1.1873^2), (0.375^2 + 1.1873^2), (0.175^2 + 0.375^2)]);
I_I2 = 1/12*m_I2 * diag([(0.3^2 + 1.1873^2), (0.3^2 + 1.1873^2), (0.3^2 + 0.3^2)]);
I_I3 = 1/12 * m_I3 * diag([(0.2^2 + 1.1873^2), (0.4^2 + 1.1873^2), (0.2^2 + 0.4^2)]);
I_I4 = 1/12*m_I4 * diag([(0.255^2 + 1.1873^2), (0.255^2 + 1.1873^2), (0.255^2 + 0.255^2)]);

%Concatenate
inertias = [I_prism I_AV I_beam1 I_beam2 I_beam3 I_I1 I_I2 I_I3 I_I4];
CoMs = [CoM_prism; CoM_AV;  CoM_beam1; CoM_beam2; CoM_beam3; CoM_I1; CoM_I2; CoM_I3; CoM_I4];
masses = [m_prism m_AV m_beam1 m_beam2 m_beam3 m_I1 m_I2 m_I3 m_I4];

%% Totals with paralell axes theorem 

% sum all of them up (except solar) with respect to CAD axes
I_tot = zeros(3,3);
for i = 1:length(masses)
    inertia_i = inertias(:, 3*(i-1)+1:3*(i-1)+3);
    CoM_i = - CoMs(i, :);
    I_tot = I_tot + inertia_i + masses(i)*(dot(CoM_i, CoM_i)*eye(3) - CoM_i' * CoM_i); %with respect to CAD
end

%Solar separately
I_tot_solar = zeros(3,3);
for i = 1:6 %add solar panels
    CoM_i = -[1.09825*cosd((i-1)*60) 1.09825*sind((i-1)*60) 0.01]; %center of mass
    I_tot_solar = I_tot_solar + I_solar + m_solar/6 *(dot(CoM_i, CoM_i)*eye(3) - CoM_i'*CoM_i); %with respect to CAD
end

%Total with respect to CAD axes
I_tot = I_tot + I_tot_solar;

%transform to body axes
CoM = - CoM;
I_tot_body = I_tot + (m_bus + m_payload) * (dot(CoM, CoM)* eye(3) - CoM' * CoM);

