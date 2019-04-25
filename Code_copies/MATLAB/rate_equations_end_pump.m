function ydot = rate_equations_end_pump(t,y,iflag,p,TT,Rp)

ydot=zeros(5,1);

t;
% USING THIS FROM NASH'S CODE DOES NOT CHANGE RESULTS
% pump_rate=interp1(TT,Rp,t,'pchip'); 

%set up the system of coupled differential equations:
ydot(1)=-y(1)/p(4)+p(2)*y(3)^2-p(3)*y(1)*y(4); 
ydot(2)=-y(2)/p(5)+y(1)/p(4); 
ydot(3)=p(1)-y(3)/p(6)+y(2)/p(5)+2*p(3)*y(1)*y(4)-2*p(2)*y(3)^2;
ydot(4)=-p(1)+((y(3)/p(6))-(p(3)*y(1))*y(4)+p(2)*y(3)^2);
ydot(5)=(-y(5)/p(14))+((p(19)*p(18)/p(11))*p(13)*((p(9)*y(3))-(p(10)*y(4)))*y(5))+((p(18)/p(11))*(y(3)/p(6))*p(15)); 