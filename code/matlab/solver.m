
Rsol=696000000;
Tsol=5777;
Rlun=1737000;
Rter=6378000;

Muter=6.7e-11*5.9737e24;
Mulun=4902000000000;

periodSol=86400*365.256;

c=299792458;

%unité astronomique
UA=149597870690;

Dter=UA;
 
alpha=0.054;

DP1=Dter*Rlun/(Rsol-Rlun);

DP3=Dter*Rlun/(Rsol*(1+alpha)-Rlun);

A1=asin(Rlun/DP1);
A2=asin(Rlun/DP3);

DP2=(DP1*tan(A1)+DP3*tan(A2))/(tan(A1)+tan(A2));
HP2=tan(A1)*(DP1-DP2);

%orbite lunaire : 
el=0.054;
wl=0;%pi/2;
Wl=0/360*2*pi;
Il=5/360*2*pi;
al=384399000;

%dephasage (pour considere les décalage initial par rapport au soleil)
fl=0;

%coeficient des fonctions :

% origine sur le point X3

p1=HP2/(DP2-DP3);
p2=HP2/(DP1-DP2);

off2=HP2+(DP2-DP3)*p2;
global ud;
global ut;
ud=al;
ut=2*pi*sqrt(al^3/Muter);

periodRatio=1;

maxT=365*2*86400;

Amax=2e-2;

global param;
param.DP3=DP3/ud;
param.DP2=DP2/ud;
param.HP2=HP2/ud;
param.DP1=DP1/ud;
param.p1=p1;
param.p2=p2;
param.off2=off2/ud;
param.Dter=Dter/ud;
param.Muter=Muter*ut*ut/ud/ud/ud;
param.Mulun=Mulun*ut*ut/ud/ud/ud;
param.al=al/ud;
param.wl=wl;
param.Wl=Wl;
param.Il=Il;
param.el=el;
param.fl=fl;
param.nl=sqrt(param.Muter/param.al^3);
param.periodRatio=periodRatio;
param.as=param.al/(param.periodRatio^(2/3));
param.ns=sqrt(param.Muter/param.as^3);
param.timeStep=300/ut;
param.prec=1e-2/ut;
param.periodSol=(periodSol*1)/ut;
param.Rter=Rter/ud;
param.Rlun=Rlun/ud;
param.maxT=maxT/ut;
param.Amax=Amax*ut*ut/ud;

reltol=1e-12;
abstol=1e-12;
optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);
tspan1 = [-1, param.maxT];
tspan = [0,-1];

X0=toCart(param.al,param.el,param.wl,param.Wl,param.Il,param.fl);
V0=getSpeed(param.al,param.el,param.wl,param.Wl,param.Il,param.fl);
x0=[X0,V0]

global sysMoon;
[x,y]=ode45(@Df, tspan, x0, optionsODE);
x0=y(end,:);
sysMoon=ode45(@Df, tspan1, x0, optionsODE);


startI=1;
endI=2;

x0=[endOBS(:,startI)',ones(1,6)]
startVal=endOBS(:,startI)'
finalVal=startOBS(:,endI)'

t0=optiNorm(3,startI)
tf=optiNorm(2,endI)

tspan=[t0,t1];

[X,Y]=ode45(@DfH,tspan,x0,optionsODE);
x0
Y(end,1:6)
finalVal

fun = @(x) payoff(t0,tf,startVal,finalVal,x);
funDelt = @(x) delta(t0,tf,startVal,finalVal,x);


disp("start")

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton',"Display","iter-detailed","MaxFunctionEvaluations",1e5);
solP=[80,-150,-10,100,6,-1]
disp("start")
solP=fminunc(fun,solP,options);
options = optimoptions('fsolve', ...
                       'Display','iter-detailed', ...
                       'FunctionTolerance',1e-8);
%solP2=fsolve(funDelt,solP,options);
solP2=solP
solP2

options = optimoptions('ga', ...
                       'Display','iter');
%solGA=ga(fun,6,options);

finalVal
value(t0,tf,startVal,finalVal,solP2)

%==========================================================================
function po=payoff(tstart,tfin,startVal,finalVal,Px)
    x0=[startVal,Px];
    tspan=[tstart,tfin];
    x0;
    reltol=1e-12;
    abstol=1e-12;
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);

    [X,Y]=ode45(@DfH,tspan,x0,optionsODE);
    Y(end,1:12);
    
    A=Y(end,7:12)/norm(Y(end,7:12));
    po=norm(finalVal-Y(end,1:6))^2;
end

%==========================================================================
function D=delta(tstart,tfin,startVal,finalVal,Px)
    x0=[startVal,Px];
    tspan=[tstart,tfin];
    x0;
    reltol=1e-12;
    abstol=1e-12;
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);

    [X,Y]=ode45(@DfH,tspan,x0,optionsODE);
    Y(end,1:12);
    
    A=Y(end,7:12)/norm(Y(end,7:12));
    D=finalVal-Y(end,1:6);
end

%==========================================================================
function X=value(tstart,tfin,startVal,finalVal,Px)
    x0=[startVal,Px];
    tspan=[tstart,tfin];
    x0;
    reltol=1e-12;
    abstol=1e-12;
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);

    [X,Y]=ode45(@DfH,tspan,x0,optionsODE);
    X=Y(end,1:6);
end

%==========================================================================
function Y=Df(t,X)
    global param;
    accL=param.Muter/norm(X(1:3))^3;
    Y=[X(4:6);-accL*X(1:3)];
end

%==========================================================================
function V=getSpeed(a,e,w,W,i,f)
   global param;
    val=sqrt(param.Muter/(a*(1-e*e)));

    CW=cos(W);
    SW=sin(W);

    CI=cos(i);

    Xnrot=sin(f+w)+e*sin(w);
    Ynrot=cos(f+w)+e*cos(w);

    x=-CW*Xnrot-SW*Ynrot*CI;
    y=-SW*Xnrot+CW*Ynrot*CI;
    z=Ynrot*sin(i);
    V=val*[x,y,z];
end

%==========================================================================
function X=toCart(a,e,w,W,i,f)
    if(e==1)
        e=0.9999999;%pour eviter les div par zéro
    end
    r=a*(1-e*e)./(1+e*cos(f));

    CT=cos(f+w);
    ST=sin(f+w);

    CW=cos(W);
    SW=sin(W);

    CI=cos(i);

    x=(CW*CT-SW*ST*CI);
    y=(SW*CT+CW*ST*CI);
    z=ST*sin(i);%on neglige la 3eme dimension
    X=r.*[x',y',z'];
end

%x : pos + vitesse + adj pos + adj speed
function dx=DfH(t,x)
    global param;
    pv=norm(x(10:12));
    u=zeros(3,1);
    if(pv>1)
        u=x(10:12)/pv*param.Amax;
    end
    accS=param.Muter/norm(x(1:3))^3;
    dx=[x(4:6);accS*(-x(1:3))+u;accS*(eye(3)-3*x(1:3)*x(1:3)'/norm(x(1:3))^2)*x(10:12);-x(7:9)];
end