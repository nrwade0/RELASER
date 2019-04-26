function ydot = TmHoYLF_rate_equations_pulsed_pump_v2_test_unit_fix_v3(t,y,iflag,p,tt,Rp)
%  y(1)=n1, y(2)=n2, y(3)=n3, y(4)=n4, y(5)=n5, y(6)=n6, y(7)=n7, y(8)=n8, y(9) = phi, y(10)=phi_out

% Based on paper "Spectroscopy and Modeling of Solid State Lanthanide
% Lasers..." (Walsh, Barnes, Petros, Yu, Singh)

% pulsed_pump_v1: Defines pump scalar pump rate, uses conditional statement in rate eq
% m-file to determine if t<=tp, in which case the pump rate is set to
% defined pump rate or 0

% pulsed_pump_v2: Defines pump pulse explicitly as a function of time, rate eq m-file
% interpolates to determine pump rate at each particular time step

% The two methods give identical results, v2 takes approx 2x longer to run

% test_unit_fix:  Fixes unit inconsistency in equations from reference
% paper (adds c*l/Lopt factor to stimulated emission terms of n7 and n8
% rate eqs, drops same factor from spontaneous emission term of photon
% density

% test_unit_fix_v2:  Same as above but uses factor of c/n instead of
% c*l/Lopt as in reference "Computational model for operation of 2um..."
% (Louchev, Urata, Saito, Wada, this paper also6 references original
% refernce paper "Spectroscopy and Modeling of Solid State Lanthanide
% Lasers...")

% test_unit_fix_v3:  Correspondances with Brian Walsh on 6/24/15 and
% 6/30/15, same as v2 except retains l/Lopt factor in second and third term
% of phi rate eq

%set up a 10 element column vector for n-dot and phi-dot:
ydot=zeros(10,1);

t;

pump_rate = interp1(tt,Rp,t,'pchip');

%set up the system of coupled differential equations:
ydot(1) = -pump_rate*(1-exp(-p(25)*p(31)*y(1))) + y(2)/p(2) + (p(16)/p(4))*y(4) + y(2)*y(8)*p(8) - y(7)*y(1)*p(9) - y(4)*y(1)*p(10) + y(2)*y(2)*p(11) + y(2)*y(7)*p(14) - y(5)*y(1)*p(15) - y(6)*y(1)*p(13) + y(3)*y(8)*p(12);
ydot(2) = -y(2)/p(2) + (p(17)/p(3))*y(3) - y(2)*y(8)*p(8) + y(7)*y(1)*p(9) + 2*y(4)*y(1)*p(10) - 2*y(2)*y(2)*p(11) - y(2)*y(7)*p(14) + y(5)*y(1)*p(15);
ydot(3) = -y(3)/p(3) + (p(18)/p(4))*y(4) + y(6)*y(1)*p(13) - y(3)*y(8)*p(12);
ydot(4) = pump_rate*(1-exp(-p(25)*p(31)*y(1))) - y(4)/p(4) - y(4)*y(1)*p(10) + y(2)*y(2)*p(11);
ydot(5) = -y(5)/p(5) + y(2)*y(7)*p(14) - y(5)*y(1)*p(15);
ydot(6) = -y(6)/p(6)    + (p(19)/p(5))*y(5) - y(6)*y(1)*p(13) + y(3)*y(8)*p(12);
ydot(7) = -y(7)/p(7) + (p(20)/p(6))*y(6) + y(2)*y(8)*p(8) - y(7)*y(1)*p(9) - y(2)*y(7)*p(14) + y(5)*y(1)*p(15) - p(23)*(1/p(33))*p(26)*(p(21)*y(7)-p(22)*y(8))*y(9);
% ydot(7) = -y(7)/p(7) + (p(20)/p(6))*y(6) + y(2)*y(8)*p(8) - y(7)*y(1)*p(9) - y(2)*y(7)*p(14) + y(5)*y(1)*p(15);
ydot(8) = y(7)/p(7) - y(2)*y(8)*p(8) + y(7)*y(1)*p(9) - y(3)*y(8)*p(12) + y(6)*y(1)*p(13) + p(23)*(1/p(33))*p(26)*(p(21)*y(7)-p(22)*y(8))*y(9);
% ydot(8) = y(7)/p(7) - y(2)*y(8)*p(8) + y(7)*y(1)*p(9) - y(3)*y(8)*p(12) + y(6)*y(1)*p(13);
ydot(9) = -y(9)/p(27) + p(23)*(1/p(33))*(p(31)/p(24))*p(26)*(p(21)*y(7)-p(22)*y(8))*y(9) + (p(31)/p(24))*(y(7)/p(7))*p(28);
ydot(10) = (log(1/p(29))*y(9))/p(30);

% ydot(1) = -Rp*(1-exp(-sigma_abs*l*y(1))) + y(2)/t2 + (B41/t4)*y(4) + y(2)*y(8)*p28 - y(7)*y(1)*p71 - y(4)*y(1)*p41 + y(2)*y(2)*p22 + y(2)*y(7)*p27 - y(5)*y(1)*p51 - y(6)*y(1)*p61 + y(3)*y(8)*p38;
% ydot(2) = -y(2)/t2 + (B32/t3)*y(3) - y(2)*y(8)*p28 + y(7)*y(1)*p71 + 2*y(4)*y(1)*p41 - 2*y(2)*y(2)*p22 - y(2)*y(7)*p27 + y(5)*y(1)*p51;
% ydot(3) = -y(3)/t3 + (B43/t4)*y(4) + y(6)*y(1)*p61 - y(3)*y(8)*p38;
% ydot(4) = Rp*(1-exp(-sigma_abs*l*y(1))) - y(4)/t4 - y(4)*y(1)*p41 + y(2)*y(2)*p22;
% ydot(5) = -y(5)/t5 + y(2)*y(7)*p27 - y(5)*y(1)*p51;
% ydot(6) = -y(6)/t6 + (B56/t5)*y(5) - y(6)*y(1)*p61 + y(3)*y(8)*p38;
% ydot(7) = -y(7)/t7 + (B67/t6)*y(6) + y(2)*y(8)*p28 - y(7)*y(1)*p71 - y(2)*y(7)*p27 + y(5)*y(1)*p51 - c*(1/n)*sigma_em*(f7*y(7)-f8*y(8))*y(9);
% ydot(8) = y(7)/t7 - y(2)*y(8)*p28 + y(7)*y(1)*p71 - y(3)*y(8)*p38 + y(6)*y(1)*p61 + c*(1/n)*sigma_em*(f7*y(7)-f8*y(8))*y(9);
% ydot(9) = -y(9)/tc + c*(1/n)*(l/Lopt)*sigma_em*(f7*y(7)-f8*y(8))*y(9) + (l/Lopt)*(y(7)/t7)*Bspon;
% ydot(10) = (log(1/Roc)*y(9))/tr;

%  all terms in rate eqs. are in [m^-3*s^-1]