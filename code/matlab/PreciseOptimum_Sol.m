

%load initial parameters :
pmeter=readtable("D:\storage\CODE\matlab\params.csv")

%bodies radius
Rsol=696000000;
Rlun=1737000;
Rter=6378000;%this value is useless for calculation
  
%gravity constant for bodies
Muter=398600441800000; %6.7e-11*5.9737e24;
Mulun=4902000000000;
MuSun=132712440018000000000;	


%unité astronomique
UA=149597870690;

%distance relative de la terre lune.
Dter=UA;
alpha=pmeter.alpha;%0.05;

%donnée sur la zone d'observation (calculer avec la lune à une distance de
%Dter du soleil
DP1=Dter*Rlun/(Rsol-Rlun);
DP3=Dter*Rlun/(Rsol*(1+alpha)-Rlun);

A1=asin(Rlun/DP1);
A2=asin(Rlun/DP3);

DP2=(DP1*tan(A1)+DP3*tan(A2))/(tan(A1)+tan(A2));
HP2=tan(A1)*(DP1-DP2);

al=384399000;

p1=HP2/(DP2-DP3);
p2=HP2/(DP1-DP2);

off2=HP2+(DP2-DP3)*p2;
%constante de normalisation.
global ud;
global ut;
ud=al;
ut=2*pi*sqrt(al^3/Muter);

%information temporelle
Synperiod=1.078;%synodic period of the Moon ( rough approximation )
periodRatio=pmeter.Dperiod;%0.92;
nbOpt=pmeter.nbOpt;%50;
maxT=365*nbOpt/20*86400;

%on normalise toutes les paramètres.
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
param.MuSun=MuSun*ut*ut/ud/ud/ud;
param.al=al/ud;
param.periodRatio=periodRatio;
param.as=param.al*(param.periodRatio^(2/3));
param.Rter=Rter/ud;
param.Rlun=Rlun/ud;
param.maxT=maxT/ut;


%load initial Condition :
dat=readtable("D:\storage\CODE\matlab\Sun_Moon_Earth_jan_2025.csv")


%on simule les trajectoires de la Terre, Lune et Soleil.
reltol=1e-12;
abstol=1e-12;
optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);
tspan1 = [-1, param.maxT];
tspan = [0,-1];

global MainSys;
moonVec=[table2array(dat(3,:)),table2array(dat(2,:))]
[x,y]=ode45(@DfTb, tspan, moonVec, optionsODE);
DfTb(0,moonVec)
x0=y(end,:);
MainSys=ode45(@DfTb, tspan1, x0, optionsODE);

t=1:0.01:param.maxT;
val=deval(MainSys,t);
X0=Origin(val);

Pe=val(7:9,:)-X0;
Pl=val(1:3,:)-X0;
hold off;
plot3(Pl(1,:),Pl(2,:),Pl(3,:))
hold on
plot3(Pe(1,:),Pe(2,:),Pe(3,:))
grid on;

xlim([-1.5,1.5])
ylim([-1.5,1.5])
zlim([-1.5,1.5])


%on opti !!!

Opt=zeros(10,nbOpt);

precOpt=zeros(5,nbOpt);

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton',"Display","off","MaxFunctionEvaluations",1e5);
x0=0.1;
fun=@(x) BestVal(x);
%[t1, y1] = fminunc(@(x) BestVal(x), x0, options);


options_fsolve = optimoptions('fsolve', ...
                            'Display','off', ...
                            'FunctionTolerance',1e-8);
%BestVal(0.31)
[t1,y1] = fsolve(@(x) BestVal(x)-100,[0],options_fsolve);
if(t1<0)
    [t1,y1] = fsolve(@(x) BestVal(x)-100,[0.5],options_fsolve);
end
[t1,y1] = fminunc(@(x) -BestVal(x), [t1], options);

%we plot the Best Val curve
N=10000;
payoff=zeros(1,(N+1));
time=(0:N)*param.maxT/N;
disp("computing trajectory")
for i=1:(N+1)
    payoff(i)=fun(time(i));
end
hold off
fig=figure()
plot(time*ut/86400,payoff);
grid on
set(gca, 'YScale', 'log')
xlim([0,maxT/86400])


%saveas(fig, "D:\storage\CODE\matlab\simpModel.png", "png");

options_fsolve = optimoptions('fsolve', ...
                            'Display','iter-detailed', ...
                            'FunctionTolerance',1e-8);
%BestVal(0.31)


T0=0;
for i=1:(N+1)
    payoff(i)=fun(time(i));
    if(payoff(i)>100 && T0==0)
        T0=time(i);
    end
end

[t1,y1] = fminunc(@(x) -BestVal(x), [T0], options);
offset=false;
if(t1<0)
    t1=t1+Synperiod;
    [t1,y1] = fminunc(@(x) -BestVal(x), [t1], options);
    offset=true;
    %[t1,y1] = fsolve(@(x) if(x>0);BestVal(x)-100;,[0.5],options_fsolve);
end

%finding dephasage:


%return;
pos=t1
%en fonction de si l'on commence par une observation descendante ou ascendante il faut calculer la première observation avant. 
if(offset)
    %opt offset
    x0=pos-5*Synperiod/6;

    [tn, yn] = fminunc(fun, x0, options);
    [out,Pm,Pz,V]=BestVal(tn);
    Opt(1,1)=tn;
    Opt(2:4,1)=Pm;
    Opt(5:7,1)=Pz;
    Opt(8:10,1)=V;
end

%on trouve les optimum dans le modèle simplifié.
for i=1:(nbOpt/2)
    %opt1
    x0=pos-Synperiod/6;

    [tn, yn] = fminunc(fun, x0, options);
    [out,Pm,Pz,V]=BestVal(tn);
    Opt(1,2*i-1+offset)=tn;
    Opt(2:4,2*i-1+offset)=Pm;
    Opt(5:7,2*i-1+offset)=Pz;
    Opt(8:10,2*i-1+offset)=V;

    %boucle if ajouté pour les même raison qu'à la ligne 174.
    if(2*i+offset<=nbOpt)
        %opt2
        x0=pos+Synperiod/6;
        [tn, yn] = fminunc(fun, x0, options);
        [out,Pm,Pz,V]=BestVal(tn);
        Opt(1,2*i+offset)=tn;
        Opt(2:4,2*i+offset)=Pm;
        Opt(5:7,2*i+offset)=Pz;
        Opt(8:10,2*i+offset)=V;
    end
    pos=pos+Synperiod;
    [pos,y]=fminunc(@(x) -BestVal(x), [pos], options);%reoptimyze to recenter the value
end

Opt

%creating result data
optiNorm=zeros(3,nbOpt);
optimum=zeros(3,nbOpt);
startOBS=zeros(6,nbOpt);
endOBS=zeros(6,nbOpt);
Xoptimum=zeros(5,nbOpt);

%essai = @(x) infoObsDE(x);

%computing real optimum from first founded approximation
for i=1:nbOpt
    i
    X0=OptiConvert(Opt(:,i));
    %solve the problem
    [X,Y]=solveProblem(X0);
    [T,TS,TF,YS,YF]=infoObsDE5(X);

    %write in data
    optiNorm(1,i)=-Y;
    optiNorm(2,i)=TS;
    optiNorm(3,i)=TF;
    optimum(1,i)=-Y*ut;
    optimum(2,i)=TS*ut;
    optimum(3,i)=TF*ut;
    
    opt=-Y*ut
    Xoptimum(:,i)=X;
    startOBS(:,i)=YS(1:6);
    endOBS(:,i)=YF(1:6);
    X0;
    X;
end
Xoptimum(:,1);
hold off;
%f=figure();
%save data
csvwrite("times.csv",optiNorm);
csvwrite("startOBS.csv",startOBS);
csvwrite("endOBS.csv",endOBS);
HSVvec=hsv(nbOpt+1);

hold off;
scatter3(0,0,0);
hold on;
for i=1:(nbOpt/2)
    plotSol(Xoptimum(:,2*i),Xoptimum(1,2*i)+1,HSVvec(2*i,:));
end
grid on;


%==========================================================================
% generate vector usable for the method ObsTimeDE from a vector computed at lines 204-207. (we fixed the speed of the device to be aligned to the
% speed of the Moon)
%==========================================================================
function X=OptiConvert(vec)
    global param;
    V=norm(vec(8:10));
    X=[vec(1),(param.DP2-param.DP3)/(param.DP1-param.DP3),sign(vec(9))*acos(vec(8)/norm(vec(8:9))),acos(vec(10)/V)];
end

%==========================================================================
% solve the full problem with an initial condition:
% first it solve with the semi major axis of the device fixed, then 
% solve the problem again using the previous solution as initial condition
% but with a fully free speed vector.
%==========================================================================
function [X,Y]=solveProblem(x0)
    disp("o");
    %disp(["pos : ",val/S*360]);
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton',"Display","off","MaxFunctionEvaluations",1e3);
    ObsTimeDE(x0);
    fun = @(x) -ObsTimeDE(x); % Negate the function to find maximum
    [X, Y] = fminunc(fun, x0, options);
    temp=makeX0DE4(X);

    x02=[X(1),X(2),temp(4:6)'];
    fun = @(x) -ObsTimeDE5(x); % Negate the function to find maximum
    [X, Y] = fminunc(fun, x02, options);
end

%==========================================================================
%compute the origin position from the position of Earth and Moon from the
%Sun
%==========================================================================
function X=Origin(vec)
    global param
    Xe=vec(7:9,:);
    Xl=vec(1:3,:);
    X=(param.Muter*Xe+param.Mulun*Xl)/(param.Muter+param.Mulun);
end

%==========================================================================
%compute the origin speed from the position of Earth and Moon from the
%Sun
%==========================================================================
function V=SpeedOrigin(vec)
    global param
    Xe=vec(10:12,:);
    Xl=vec(4:6,:);
    V=(param.Muter*Xe+param.Mulun*Xl)/(param.Muter+param.Mulun);
end

%==========================================================================
% dynamics of the Earth Moon Sun system (the sun is fixed at the origin)
%==========================================================================
function Y=DfTb(t,X)
    global param;
    %global SunSys;

    %Si=deval(SunSys,t);
    %Sp=Si(1:3);
    %acc0S=param.MuSun/norm(Sp)^3;
    accLS=param.MuSun/norm(X(1:3))^3;
    accES=param.MuSun/norm(X(7:9))^3;

    accL=param.Muter/norm(X(1:3)-X(7:9))^3;
    accT=param.Mulun/norm(X(7:9)-X(1:3))^3;

    Y=[X(4:6);-accL*(X(1:3)-X(7:9))-accLS*X(1:3);X(10:12);-accT*(X(7:9)-X(1:3))-accES*X(7:9)];
end

%==========================================================================
%compute the positions of Earth and Moon in the coordinate system using the
%solved system
%==========================================================================
function [X,V,Xe,Ve]=EarthMoonPos(t)
    global MainSys;
    Y=deval(MainSys,t);
    X0=Origin(Y);
    V0=SpeedOrigin(Y);
    X=Y(1:3,:)-X0;
    V=Y(4:6,:)-V0;
    Xe=Y(7:9,:)-X0;
    Ve=Y(10:12,:)-V0;
end

%==========================================================================
%position of the Sun relative to the origin (as the sun is fixed at the
% the origin of the Earth Moon sun System, the Sun position is simply the
% negative of the origin position)
%==========================================================================
function Xst=SolPos(t)
    global MainSys;
    
    Y=deval(MainSys,t);
    X0=Origin(Y);
    Xst=-X0;
end

%==========================================================================
%position of the Sun relative to the origin (as the sun is fixed at the
% the origin of the Earth Moon sun System, the Sun position is simply the
% negative of the origin position)
%==========================================================================
function Vst=SolSpeed(t)
    global MainSys;
    Y=deval(MainSys,t);
    V0=SpeedOrigin(Y);
    Vst=-V0;
end

%==========================================================================
% objective function in the simplified model (use the relative velocity of
% the observation zone compared to the speed the device would have at this
% altitude.
%==========================================================================
function [out,Pm,Pz,V]=BestVal(t)
    global param;
    %pos du soleil
    
    Xst=SolPos(t);
    Vst=SolSpeed(t);

    %pos de la lune
    [posLun,Vl,posEarth,Ve]=EarthMoonPos(t);
    %step 2: calculer la position de la zone d'observation (facile)
    P2=posLun+param.DP2/param.Dter*(posLun-Xst);

    %step 3: calculer la vitesse de la zone d'observation (facile)
    VP2=Vl+param.DP2/param.Dter*(Vl-Vst);
    

    %step 4: calculer la vitesse du satelite:
    VsN=abs(sqrt(param.Muter*(2/norm(P2)-1/param.as)));

    %retourner le payoff (très facile)
    %out=(norm(posLun)-norm(P2))^2;
    out=(norm(VP2)-VsN)^2;
    Test=param.Muter*(2/norm(P2)-1/param.as);
    Pm=posLun;
    Pz=P2;
    V=VP2;%getSpeed(param.al,param.el,param.wl,param.Wl,param.Il,f);
end

%==========================================================================
% this method check if the device is inside the observation zone 
%==========================================================================
function [value, isterminal, direction]=inObsDE(t,X)
    global param;
    global sysMoonEarth;
    [Xl,Vl,Xe,Ve]=EarthMoonPos(t);
    Ps=SolPos(t);
    %Xl=X(7:9)';
    P3=Xl+(param.DP3/param.Dter)*(Xl-Ps);%*sqrt(param.Dter*(param.Dter+2*Xl(1))))*[1,Xl(2)/(param.Dter+Xl(1)),0];%[1,X(2)/(X(1)+param.Dter),0];
    P1=Xl+(param.DP1/param.Dter)*(Xl-Ps);
    

    Dir=P1-P3;

    PX=dot(X(1:3)'-P3',Dir)/norm(Dir)^2*Dir;
    PY=X(1:3)-P3-PX;
    %POS=[X(1)-X(7),X(2)-X(8),X(3)-X(9)]-P3;
    d=norm(PY);%sqrt(POS(2)*POS(2)+POS(3)*POS(3));%norm(POS());
    sx=norm(PX);%POS(1);
    sx=sx*(param.DP1-param.DP3)/norm(Dir);
    out=d<param.p1*sx && -d<param.p1*sx && d<(param.off2-param.p2*sx) && -d<(param.off2-param.p2*sx);

    value=out-0.5;
    isterminal=1;
    direction = -1;
end

%==========================================================================
% this method compute the dynamics of the deivce among the Moon and Earth,
% the Sun gravity variation is neglected as the device stay inside a sphere
% of two Earth Moon distance from the Earth.
%==========================================================================
function Y=Dfl(t,X)
    global param;
    global sysMoonEarth;

    [Xl,Vl,Xe,Ve]=EarthMoonPos(t);

    accS=param.Muter/norm(X(1:3)-Xe)^3;
    
    %gravité de la lune sur l'objet
    accSl=param.Mulun/norm(X(1:3)-Xl)^3;

    %accL=param.Muter/norm(X(7:9))^3;
    Y=[X(4:6);-accS*(X(1:3)-Xe)-accSl*(X(1:3)-Xl)];
end

%==========================================================================
% plot the trajectory of the device during the observation
%==========================================================================
function T=plotSol(X,tmax,col)
    global param;
    %after :
    x0=makeX0DE5(X);%;[R,V,Xl,Vl]';
    
    if(inObsDE(X(1),x0)<0)
        disp("probleme");
        T=0;
        return;
    end
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    
    %optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol,'Events', @inObsDE);
    tspan = [X(1), tmax];
    [t, Y] = ode45(@Dfl, tspan, x0, optionsODE);

    plot3(Y(:,1),Y(:,2),Y(:,3),'Color', col);

    tspan = [X(1), -1];
    [t, Y] = ode45(@Dfl, tspan, x0, optionsODE);

    plot3(Y(:,1),Y(:,2),Y(:,3),'Color', col);

end

%==========================================================================
%compute the observation time for the given vector, the vector consist of 5
%value : 
% the time when the device cross the center of the observation zone.
% the cross point of the device on the center line of the observation zone.
% two angles corresponding to the direction of the speed:
% the speed magnitude
%==========================================================================
function T=ObsTimeDE5(X)
    global param;
    %after :

    x0=makeX0DE5(X);
    if(inObsDE(X(1),x0)<0)
        disp("probleme");
        T=0;
        return;
    end
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol,'Events', @inObsDE);
    tspan = [X(1), X(1)+1];
    [t, Y,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=TE;
    tspan = [X(1), X(1)-1];
    [t, Y,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=T-TE;

end

%==========================================================================
%compute the observation time for the given vector, the vector consist of 5
%value : 
% the time when the device cross the center of the observation zone.
% the cross point of the device on the center line of the observation zone.
% two angles corresponding to the direction of the speed:
% the speed is deteminded witht the parameter as.
%==========================================================================
function T=ObsTimeDE(X)
    global param;
    %after :

    x0=makeX0DE4(X);
    if(inObsDE(X(1),x0)<0)
        disp("probleme");
        T=0;
        return;
    end
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol,'Events', @inObsDE);
    tspan = [X(1), X(1)+1];
    [t, Y,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=TE;
    tspan = [X(1), X(1)-1];
    [t, Y,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=T-TE;

end

%==========================================================================
% generate the initial position and velocity of the device from the vector
% given in input of the method ObsTimeDE (it is separated from this method
% because it could be used in multiple place)
%==========================================================================
function x0=makeX0DE4(X)
    global param;
    [Xl,Vl,Xe,Ve]=EarthMoonPos(X(1));%toCart(param.al,param.el,param.wl,param.Wl,param.Il,X(1));
    %XP3=[param.DP3,0,0];
    %XP1=[param.DP1,0,0];
    Xts=SolPos(X(1));
    XP3=Xl'+(param.DP3/param.Dter)*(Xl'-Xts');
    XP1=Xl'+(param.DP1/param.Dter)*(Xl'-Xts');

    %radius:
    R=(1-X(2))*(XP3)+X(2)*XP1;
    Rn=norm(R);
    Vn=1*sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*sin(X(4)),sin(X(4))*sin(X(3)),cos(X(4))]*Vn;

    x0=[R,V]';
end

%==========================================================================
% generate the initial position and velocity of the device from the vector
% given in input of the method ObsTimeDE5 (it is separated from this method
% because it could be used in multiple place)
%==========================================================================
function x0=makeX0DE5(X)
    global param;
    
    [Xl,Vl,Xe,Ve]=EarthMoonPos(X(1));%toCart(param.al,param.el,param.wl,param.Wl,param.Il,X(1));
    %XP3=[param.DP3,0,0];
    %XP1=[param.DP1,0,0];
    Xts=SolPos(X(1));
    XP3=Xl'+(param.DP3/param.Dter)*(Xl'-Xts');
    XP1=Xl'+(param.DP1/param.Dter)*(Xl'-Xts');

    %radius:
    R=(1-X(2))*(XP3)+X(2)*XP1;
    Rn=norm(R);
    %Vn=X(5);%1*sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[X(3),X(4),X(5)];%[cos(X(3))*sin(X(4)),sin(X(4))*sin(X(3)),cos(X(4))]*Vn;

    x0=[R,V]';
end

%==========================================================================
% do pretty much the same things as the method ObsTimeDE5 except it return
% more informations about the observations (used after the optimisation to
% fill the csv)
%==========================================================================
function [T,TS,TF,YS,YF,V]=infoObsDE5(X)
    global param;
    global sysMoonEarth;
    %after :

    x0=makeX0DE5(X);
    if(inObsDE(X(1),x0)<0)
        disp("probleme");
        T=0;
        return;
    end
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol,'Events', @inObsDE);
    tspan = [X(1), X(1)+1];
    [t, Y,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=TE;
    TF=TE;
    YF=YE;
    tspan = [X(1), X(1)-1];
    [t, Y,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=T-TE;
    TS=TE;
    YS=YE;
end
