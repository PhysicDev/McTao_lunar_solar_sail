
Rsol=696000000;
Tsol=5777;
Rlun=1737000;
Rter=6500000;

Muter=6.7e-11*5.9737e24;

c=299792458;

%unité astronomique
UA=149597870690;

Dter=UA;

alpha=0.1;

DP1=Dter*Rlun/(Rsol-Rlun);

DP3=Dter*Rlun/(Rsol*(1+alpha)-Rlun);

A1=asin(Rlun/DP1);
A2=asin(Rlun/DP3);

DP2=(DP1*tan(A1)+DP3*tan(A2))/(tan(A1)+tan(A2));
HP2=tan(A1)*(DP1-DP2);

%orbite lunaire : 
el=0.0+1e-8;%54;
wl=0+1e-8;%pi/2;
al=384399000;
Wl=0+1e-8;
il=pi/4+1e-8;
%coeficient des fonctions :

% origine sur le point X3

p1=HP2/(DP2-DP3);
p2=HP2/(DP1-DP2);

off2=HP2+(DP2-DP3)*p2;

ud=al;
ut=2*pi*sqrt(al^3/Muter);

periodRatio=0.5;

global param;
param.DP3=DP3/ud;
param.DP1=DP1/ud;
param.p1=p1;
param.p2=p2;
param.off2=off2/ud;
param.Muter=Muter*ut*ut/ud/ud/ud;
param.al=al/ud;
param.wl=wl;
param.el=el;
param.Wl=Wl;
param.il=il;
param.nl=sqrt(param.Muter/param.al^3);
param.periodRatio=periodRatio;
param.as=param.al/(param.periodRatio^(2/3));
param.ns=sqrt(param.Muter/param.as^3);
param.timeStep=300/ut;
param.prec=1e-2/ut;


Xtest=[4*pi/6,0.5,150/360*2*pi,0];
T=ObsTimeDE(Xtest)*ut





x0=[2*pi/3,0.5,210/360*2*pi,0];
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton',"Display","iter-detailed","MaxFunctionEvaluations",1500);
fun = @(x) -ObsTimeDE(x); % Negate the function to find maximum
[x, y] = fminunc(fun, x0, options);
x0
disp(x);
disp(-y*ut);
disp(ObsTimeDE(x)*ut);



function out=InObs(x)
    global param;
    d=sqrt(x(2)*x(2)+x(3)*x(3));
    out=(d<param.p1*x(1) & -d<param.p1*x(1) & d<(param.off2-param.p2*x(1)) & -d<(param.off2-param.p2*x(1)));
end

function [value, isterminal, direction]=inObsDE(t,X)
    global param;
    value=InObs([X(1)-X(7)-param.DP3,X(2)-X(8),X(3)-X(9)])-0.5;
    isterminal=1;
    direction = -1;
end

function Y=Df(t,X)
    global param;
    accS=param.Muter/norm(X(1:3))^3;
    accL=param.Muter/norm(X(7:9))^3;
    Y=[X(4:6);-accS*X(1:3);X(10:12);-accL*X(7:9)];
end


function out=plotOrbit(X)
    global param;
    %after :
    while(X(1)<=0)
        X(1)=X(1)+2*pi
    end
    [Xl,Vl]=keplerian2ijk(param.al,param.el,param.il,param.Wl,param.wl,X(1));%toCart(param.al,param.el,param.wl,X(1));
    XP3=[param.DP3,0,0];
    XP1=[param.DP1,0,0];
    
    %radius:
    R=Xl'+XP3*(1-X(2))+XP1*X(2);
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*cos(X(4)),cos(X(4))*sin(X(3)),sin(X(4))]*Vn;
    
    x0=[R,V,Xl',Vl']';
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);
    tspan = [0, 1];
    [t, Y] = ode45(@Df, tspan, x0, optionsODE);
    plot(Y(:,1),Y(:,2));
    grid on;
    hold off;
end

function T=ObsTimeDE(X)
    global param;
    %after :
        


    while(X(1)<=0)
        X(1)=X(1)+2*pi
    end
    [Xl,Vl]=keplerian2ijk(param.al,param.el,param.il/pi*180,param.Wl/pi*180,param.wl/pi*180,X(1)/pi*180);%toCart(param.al,param.el,param.wl,X(1));
    XP3=[param.DP3,0,0];
    XP1=[param.DP1,0,0];
    %update of Vl because didn't find a way to indicate Mu.
    %radius:
    R=Xl'+XP3*(1-X(2))+XP1*X(2);
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*cos(X(4)),cos(X(4))*sin(X(3)),sin(X(4))]*Vn;
    Vl=Vl/norm(Vl)*sqrt(param.Muter*(2/norm(Xl)-1/param.al));
    x0=[R,V,Xl',Vl']';
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol,'Events', @inObsDE);
    tspan = [0, 0.5];
    [t, Y,TE,YE,IE] = ode45(@Df, tspan, x0, optionsODE);
    T=TE;
    tspan = [0, -0.5];
    [t, Y,TE,YE,IE] = ode45(@Df, tspan, x0, optionsODE);
    T=T-TE;

end

function f=posInOrb(e,n,t)
    M=t*n;
    options = optimoptions('fsolve', ...
                            'Display','off', ...
                            'FunctionTolerance',1e-8);
    E=fsolve(@(x) (x-e*sin(x)-M),M,options);
    f=2*atan(sqrt((1+e)/(1-e))*tan(E/2));%sqrt((1+e)/(1-e))
end

function T=perigeeTime(e,n,f)
    E=2*atan(sqrt((1-e)/(1+e))*tan(f/2));
    M=E-e*sin(E);
    T=M/n;
end