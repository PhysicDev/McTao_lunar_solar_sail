
Rsol=696000000;
Tsol=5777;
Rlun=1737000;
Rter=6378000;
  
Muter=398600441800000; %6.7e-11*5.9737e24;
Mulun=4902000000000;
MuSun=132712440018000000000;	

periodSol=86400*365.256;

c=299792458;

%unité astronomique
UA=149597870690;

Dter=UA;
 
alpha=0.005;

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


Synperiod=1.078;%synodic period of the Moon ( rough approximation )
periodRatio=0.92;
nbOpt=50;
maxT=365*nbOpt/20*86400;

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
param.wl=wl;
param.Wl=Wl;
param.Il=Il;
param.el=el;
param.fl=fl;
param.nl=sqrt(param.Muter/param.al^3);
param.periodRatio=periodRatio;
param.as=param.al*(param.periodRatio^(2/3));
param.ns=sqrt(param.Muter/param.as^3);
param.timeStep=300/ut;
param.prec=1e-2/ut;
param.periodSol=(periodSol*1)/ut;
param.Rter=Rter/ud;
param.Rlun=Rlun/ud;
param.maxT=maxT/ut;


%load initial Condition :
dat=readtable("D:\storage\CODE\matlab\Sun_Moon_Earth_jan_2025.csv")



reltol=1e-12;
abstol=1e-12;
optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);
tspan1 = [-1, param.maxT];
tspan = [0,-1];


%global SunSys

%simulate Sun movement : 
%sunVec=[table2array(dat(1,:))]
%[x,y]=ode45(@DfSun, tspan, sunVec, optionsODE);
%x0=y(end,:);
%SunSys=ode45(@DfSun, tspan1, x0, optionsODE);

%computing moon and earth movement : 

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
for i=1:(nbOpt/2)
    %opt1
    x0=pos-Synperiod/6;

    [tn, yn] = fminunc(fun, x0, options);
    [out,Pm,Pz,V]=BestVal(tn);
    Opt(1,2*i-1+offset)=tn;
    Opt(2:4,2*i-1+offset)=Pm;
    Opt(5:7,2*i-1+offset)=Pz;
    Opt(8:10,2*i-1+offset)=V;

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

optiNorm=zeros(3,nbOpt);
optimum=zeros(3,nbOpt);
startOBS=zeros(6,nbOpt);
endOBS=zeros(6,nbOpt);
Xoptimum=zeros(5,nbOpt);

%essai = @(x) infoObsDE(x);

for i=1:nbOpt
    i
    X0=OptiConvert(Opt(:,i));
    [X,Y]=solveProblem(X0);
    [T,TS,TF,YS,YF]=infoObsDE5(X);

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
csvwrite("times.csv",optiNorm);
csvwrite("startOBS.csv",startOBS);
csvwrite("endOBS.csv",endOBS);
HSVvec=hsv(nbOpt+1);

%zoneMove(Xoptimum(:,1),150);
%hold on;
%for i=1:nbOpt
%    if(mod(i,2)==1)
%        visu3D(Xoptimum(:,i),HSVvec(i,:));
%    end
%end
%f=figure();
%hold on;

%for i=1:nbOpt
%    plotSol(Xoptimum(:,i),optimum(1,end)+0.5);
%end
%grid on;
%legend(1:nbOpt)

function out=zoneMove(X,N)
    global param;
    global ut;
    %after :

    [Xl,Vl]=MoonPos(X(1));%toCart(param.al,param.el,param.wl,param.Wl,param.Il,X(1));
    %XP3=[param.DP3,0,0];
    %XP1=[param.DP1,0,0];
    Xts=SolPos(X(1));
    XP3=(param.DP3/param.Dter)*([Xl(1),Xl(2),Xl(3)]-Xts);
    XP1=(param.DP1/param.Dter)*([Xl(1),Xl(2),Xl(3)]-Xts);

    %radius:
    R=Xl+(1-X(2))*(XP3)+X(2)*XP1;
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*sin(X(4)),sin(X(4))*sin(X(3)),cos(X(4))]*Vn;

    x0=[R,V,Xl,Vl]';
    
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
    TS=TE;
    tspan = [X(1), X(1)-1];
    [t, Y,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=T-TE;


    optionsODE = odeset('RelTol', reltol, 'AbsTol',abstol);
    X0=YE;
    tspan = [TE,TS];
    Sol = ode45(@Dfl, tspan, X0, optionsODE);
    tspan = [X(1),TS];
    Sol2 = ode45(@Dfl, tspan, x0, optionsODE);

    timeVec=(0:(N-1))/(N-1)*(TE-TS)+TS;
    Yint=deval(Sol,timeVec);
    Psol=SolPos(timeVec');
    Rl=[Yint(7,:)',Yint(8,:)',Yint(9,:)'];
    for i=1:(N-1)
        disp(i);
        f=figure("Visible","off");
        drawObs(Rl(i,:)',Psol(i,:));
        hold on;
        scatter([Yint(1,i)],[Yint(2,i)]);
        plot(Yint(1,:),Yint(2,:));
        title("time "+(timeVec(i)-TE)*ut/3600+" h");
        hold off;

        saveas(f, "D:\storage\CODE\matlab\framesXY\frame_"+(10000+i), "png");
        close(f);
        f=figure("Visible","off");
        drawObsXZ(Rl(i,:)',Psol(i,:));
        hold on;
        scatter([Yint(1,i)],[Yint(3,i)]);
        plot(Yint(1,:),Yint(3,:));
        title("time "+(timeVec(i)-TE)*ut/3600+" h");
        hold off;

        saveas(f, "D:\storage\CODE\matlab\framesXZ\frame_"+(10000+i), "png");
        close(f);
    end
end

%==========================================================================
function out=visu3D(X,col)
    global param;
    %after :

    [Xl,Vl]=MoonPos(X(1));%toCart(param.al,param.el,param.wl,param.Wl,param.Il,X(1));
    %XP3=[param.DP3,0,0];
    %XP1=[param.DP1,0,0];
    Xts=SolPos(X(1));
    XP3=(param.DP3/param.Dter)*([Xl(1),Xl(2),Xl(3)]-Xts);
    XP1=(param.DP1/param.Dter)*([Xl(1),Xl(2),Xl(3)]-Xts);

    %radius:
    R=Xl+(1-X(2))*(XP3)+X(2)*XP1;
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*sin(X(4)),sin(X(4))*sin(X(3)),cos(X(4))]*Vn;

    x0=[R,V,Xl,Vl]';
    
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
    TS=TE;
    tspan = [X(1), X(1)-1];
    [t, Y,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=T-TE;


    optionsODE = odeset('RelTol', reltol, 'AbsTol',abstol);
    X0=YE;
    tspan = [TE,TS];
    Sol = ode45(@Dfl, tspan, X0, optionsODE);
    tspan = [X(1),TS];
    Sol2 = ode45(@Dfl, tspan, x0, optionsODE);

    %print moon orbit;
    %disp("okok");
    %hold off;
    plotMoon3D(1);
    hold on;
    scatter3([0],[0],[0]);
    %scatter3([],[],[]);
    n=100;

    timeVec=(0:(n-1))/(n-1)*(TE-TS)+TS;
    Yint=deval(Sol,timeVec);
    Psol=SolPos(timeVec');
    Rl=[Yint(7,:)',Yint(8,:)',Yint(9,:)'];
    %size([Yint(7,:)',Yint(8,:)',Yint(9,:)'])
    %size(Yint)
    %size(Psol)
    P3=Rl+(param.DP3/param.Dter)*(Rl-Psol);
    P1=Rl+(param.DP1/param.Dter)*(Rl-Psol);
    P2=Rl+(param.DP2/param.Dter)*(Rl-Psol);
    %disp("are you ploting ?")
    %ax = gca;
    %ax.ColorOrder = HSVvec;
    %ax.ColorOrderIndex=1;
    plot3([P1(:,1),P3(:,1)]',[P1(:,2),P3(:,2)]',[P1(:,3),P3(:,3)]',"Color","black");

    plot3(Yint(1,:)',Yint(2,:)',Yint(3,:)',"-o","Color",col,'MarkerSize',5);
    plot3(Yint(7,:)',Yint(8,:)',Yint(9,:)',"-x","Color",col,'MarkerSize',5);
    %plot3(Y(1,:)',Y(2,:)',Y(3,:)',"Color","blue")
    %hold off;
end

%==========================================================================
function X=OptiConvert(vec)
    global param;
    V=norm(vec(8:10));
    X=[vec(1),(param.DP2-param.DP3)/(param.DP1-param.DP3),sign(vec(9))*acos(vec(8)/norm(vec(8:9))),acos(vec(10)/V)];
end

%==========================================================================
function [X,Y]=solveProblem(x0)
    disp("o");
    %disp(["pos : ",val/S*360]);
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton',"Display","off","MaxFunctionEvaluations",1e3);
    ObsTimeDE(x0)
    fun = @(x) -ObsTimeDE(x); % Negate the function to find maximum
    [X, Y] = fminunc(fun, x0, options);
    temp=makeX0DE4(X);
    x02=[X,norm(temp(4:6))]
    fun = @(x) -ObsTimeDE5(x); % Negate the function to find maximum
    [X, Y] = fminunc(fun, x02, options);
end

%==========================================================================
function dvt=deltaVtime(t1,t2,sol1,sol2,Stime,Etime)
    global param;
    t1;
    t2;
    Stime;
    Etime;
    if((t2-Etime)>0)
        dvt=1000;
        return;
    end
    if(t1-Stime<0)
        dvt=2000;
        return;
    end
    if(t1>t2)
        dvt=3000;
        return;
    end

    Ys1=deval(sol1,t1-Stime)';
    Ys2=deval(sol2,t2-Etime)';
    
    Ps1=Ys1(:,1:3);
    Vs1=Ys1(:,4:6);
           
    Ps2=Ys2(:,1:3);
    Vs2=Ys2(:,4:6);


    dt=(t2-t1)';

    [Vi,Vf]=solvelambert(Ps1,Ps2,dt,ones(length(dt), 1),param.Muter);
    dvt=vecnorm((Vi-Vs1)')+vecnorm((Vf-Vs2)');
end

%==========================================================================
function T=perigeeTime(e,n,f)
    E=2*atan(sqrt((1-e)/(1+e))*tan(f/2));
    M=E-e*sin(E);
    T=M/n;
end

%==========================================================================
function f=posInOrb(e,n,t)
    M=t*n;
    options = optimoptions('fsolve', ...
                            'Display','off', ...
                            'FunctionTolerance',1e-8);
    E=fsolve(@(x) (x-e*sin(x)-M),M,options);
    f=2*atan(sqrt((1+e)/(1-e))*tan(E/2));%sqrt((1+e)/(1-e))
end

%==========================================================================
function Y=DfSun(t,X)%compute dynamics of the Sun (ie movement of the sun relative to the Earth-Moon barycenter)
    global param;
    accSys=param.MuSun/norm(X(1:3))^3;
    Y=[X(4:6);-accSys*X(1:3)];
end

%==========================================================================
function Y=Df(t,X)
    global param;
    accL=param.Muter/norm(X(1:3))^3;
    Y=[X(4:6);-accL*X(1:3)];
end

%==========================================================================
function X=Origin(vec)
    global param
    Xe=vec(7:9,:);
    Xl=vec(1:3,:);
    X=(param.Muter*Xe+param.Mulun*Xl)/(param.Muter+param.Mulun);
end

%==========================================================================
function V=SpeedOrigin(vec)
    global param
    Xe=vec(10:12,:);
    Xl=vec(4:6,:);
    V=(param.Muter*Xe+param.Mulun*Xl)/(param.Muter+param.Mulun);
end

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

%==========================================================================
function drawObsXZ(Xl,Xst)
    %on desine en 2D donc on neglige l'axe z pour l'instant
    global param;
    Xls=Xl'-Xst;%[param.Dter+Xl(1),Xl(2),Xl(3)];
    XP3=(param.DP3/param.Dter)*Xls;
    XP1=(param.DP1/param.Dter)*Xls;
    XP2x=(param.DP2/param.Dter)*Xls;
    XP2y=(XP2x/norm(XP2x))*[0,0,-1;0,1,0;1,0,0]*param.HP2;
    poly=[XP1+Xl';Xl'+XP2x+XP2y;Xl'+XP3;Xl'+XP2x-XP2y;XP1+Xl'];
    plot(poly(:,1),poly(:,3));


    xlim([min(poly(:,1)),max(poly(:,1))])
    ylim([min(poly(:,3)),max(poly(:,3))])


end

%==========================================================================
function drawObs(Xl,Xst)
    %on desine en 2D donc on neglige l'axe z pour l'instant
    global param;
    Xls=Xl'-Xst;%[param.Dter+Xl(1),Xl(2),Xl(3)];
    XP3=(param.DP3/param.Dter)*Xls;
    XP1=(param.DP1/param.Dter)*Xls;
    XP2x=(param.DP2/param.Dter)*Xls;
    XP2y=(XP2x/norm(XP2x))*[0,-1,0;1,0,0;0,0,1]*param.HP2;
    poly=[XP1+Xl';Xl'+XP2x+XP2y;Xl'+XP3;Xl'+XP2x-XP2y;XP1+Xl'];
    plot(poly(:,1),poly(:,2));

    xlim([min(poly(:,1)),max(poly(:,1))])
    ylim([min(poly(:,2)),max(poly(:,2))])
end

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
function Xst=SolPos(t)
    global MainSys;
    
    Y=deval(MainSys,t);
    X0=Origin(Y);
    Xst=-X0;
end

%==========================================================================
function Vst=SolSpeed(t)
    global MainSys;
    Y=deval(MainSys,t);
    V0=SpeedOrigin(Y);
    Vst=-V0;
end

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
function plotMoon(t)
    global param;
    Xl=toCart(param.al,param.el,param.wl,param.Wl,param.Il,param.fl);
    Vl=getSpeed(param.al,param.el,param.wl,param.Wl,param.Il,param.fl);
    x0=[Xl,Vl];
    tspan=[0,1];
    [time, Y] = ode45(@Df, tspan, x0);
    hold on;
    plot( Y(:,1), Y(:,2));
    grid on;
    grid minor;

    T=perigeeTime(param.el,param.nl,param.fl);
    f=posInOrb(param.el,param.nl,t-T);
    posLun=toCart(param.al,param.el,param.wl,param.Wl,param.Il,f);

    timeV=(0:100)/100*2*pi;
    plot(cos(timeV)*param.Rter,sin(timeV)*param.Rter)
    plot(cos(timeV)*param.Rter+posLun(1),sin(timeV)*param.Rter+posLun(2))
end 

%==========================================================================
function plotMoon3D(t)
    global param;
    Xl=toCart(param.al,param.el,param.wl,param.Wl,param.Il,param.fl);
    Vl=getSpeed(param.al,param.el,param.wl,param.Wl,param.Il,param.fl);
    x0=[Xl,Vl];
    tspan=[0,t];
    [time, Y] = ode45(@Df, tspan, x0);
    plot3( Y(:,1), Y(:,2),Y(:,3));
    grid on;
    grid minor;
end 

%==========================================================================
function plotOrbitXY(a,e,w,W,i)
    angle=-pi:0.01:pi;
    X=toCart(a,e,w,W,i,angle);
    plot(X(:,1),X(:,2));
end

%==========================================================================
function [a,e,w,W,i,f] = getOrbital(r, v, mu)

    % Specific angular momentum
    h = cross(r, v);
    h_norm = norm(h);
    
    % Inclination
    i = acos(h(3) / h_norm);
    
    % Node line
    K = [0, 0, 1];
    N = cross(K, h);
    N_norm = norm(N);
    
    % Right Ascension of the Ascending Node (RAAN)
    if N_norm ~= 0
        RAAN = acos(N(1) / N_norm);
        if N(2) < 0
            RAAN = 2 * pi - RAAN;
        end
    else
        RAAN = 0;
    end
    
    % Eccentricity vector
    v_r = dot(r, v) / norm(r);
    e = (1 / mu) * ((norm(v)^2 - mu / norm(r)) * r - norm(r) * v_r * v);
    e_norm = norm(e);
    
    % Argument of Periapsis
    if N_norm ~= 0
        if e_norm > 1e-10 % Threshold to avoid division by zero
            omega = acos(dot(N, e) / (N_norm * e_norm));
            if e(3) < 0
                omega = 2 * pi - omega;
            end
        else
            omega = 0;
        end
    else
        omega = 0;
    end
    
    % True Anomaly
    if e_norm > 1e-10 % Threshold to avoid division by zero
        f = acos(dot(e, r) / (e_norm * norm(r)));
        if dot(r, v) < 0
            f = 2 * pi - f;
        end
    else
        cp = cross(N, r);
        if cp(3) >= 0
            f = acos(dot(N, r) / (N_norm * norm(r)));
        else
            f = 2 * pi - acos(dot(N, r) / (N_norm * norm(r)));
        end
    end
    
    % Semi-major axis
    a = 1 / ((2 / norm(r)) - (norm(v)^2 / mu));
    
    % Output structure
    e=e_norm;
    w=omega;
    W=RAAN;
    i=i;
    f=f;
end

%==========================================================================
function P3=getXP3(Xl)
    global param;
    P3=(param.DP3/param.Dter)*[param.Dter+Xl(1),Xl(2),Xl(3)];%*sqrt(param.Dter*(param.Dter+2*Xl(1))))*[1,Xl(2)/(param.Dter+Xl(1)),0];%[1,X(2)/(X(1)+param.Dter),0];
end

%==========================================================================
function out=InObs(x)
    global param;
    d=sqrt(x(2)*x(2)+x(3)*x(3));
    out=all([d<param.p1*x(1) && -d<param.p1*x(1) && d<(param.off2-param.p2*x(1)) && -d<(param.off2-param.p2*x(1))]);
end

%==========================================================================
function plot_inter_t(t,y,step)
    global ut;
    global param;
    frames=(t(end)-t(1))/step;
    disp(["frames  :  ",frames]);%length(DifY(:,1))])

    Ct=t(1);
    lvec=1:length(t)

    for i=1:frames%length(DifY(:,1))
        
        pos=lvec(t<=Ct);
        pos=pos(end);
        sous_t=t(pos);
        next_t=t(pos+1);

        Xs=SolPos(Ct);

        %on crée
        lambda=(Ct-next_t)/(sous_t-next_t);

        Y=y(pos,:)*lambda+y(pos+1,:)*(1-lambda);

        disp(["progress",i/frames]);
        fig=figure('Visible', 'off');
        hold on;
        %plot(Y(:,1),Y(:,2));
        %plot(Y(:,7),Y(:,8));
        %scatter([DifY(i,1),DifY(i,7)],[DifY(i,2),DifY(i,8)])
        scatter([Y(1)],[Y(2)])
        drawObs(Y(7:9)',Xs);
        grid on;

        title("T="+Ct*ut/3600+"h");
        hold off;
        saveas(fig, "D:\storage\CODE\matlab\framesXY\frame_"+i+".png","png");
        close(fig);

        fig=figure('Visible', 'off');
        hold on;
        %plot(Y(:,1),Y(:,3));
        %plot(Y(:,7),Y(:,9));
        %scatter([DifY(i,1),DifY(i,7)],[DifY(i,2),DifY(i,8)])
        scatter([Y(1)],[Y(3)])
        drawObsXZ(Y(7:9)',Xs);
        grid on;

        title("T="+Ct*ut/3600+"h");
        hold off;
        saveas(fig, "D:\storage\CODE\matlab\framesXZ\frame_"+i+".png","png");
        close(fig);

        Ct=Ct+step;
    end
end

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
function out=plotOrbit(X)
    global param;
    global ut;
    %after :

    [Xl,Vl]=MoonPos(X(1));%toCart(param.al,param.el,param.wl,param.Wl,param.Il,X(1));
    %XP3=[param.DP3,0,0];
    %XP1=[param.DP1,0,0];
    Xts=SolPos(X(1));
    XP3=(param.DP3/param.Dter)*([Xl(1),Xl(2),Xl(3)]-Xts);
    XP1=(param.DP1/param.Dter)*([Xl(1),Xl(2),Xl(3)]-Xts);

    %radius:
    R=Xl+(1-X(2))*(XP3)+X(2)*XP1;
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*sin(X(4)),sin(X(4))*sin(X(3)),cos(X(4))]*Vn;

    x0=[R,V,Xl,Vl]';
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    %fig=figure();
    fig=figure('Visible', 'off');
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);
    tspan = [0, 2];
    [t, Y] = ode45(@Dfl, tspan, x0, optionsODE);
    hold on;
    plot(Y(:,1),Y(:,2));
    plot(Y(:,7),Y(:,8));
    xlim([-1,2]);
    ylim([-1,1]);
    grid on;
    scatter([x0(1),x0(7)],[x0(2),x0(8)])
    %drawObs(x0(7:9));
    hold off;
    saveas(fig, "D:\storage\CODE\matlab\framesXY\orbitXY_"+param.Wl+".png","png");
    fig=figure('Visible', 'off');
    %fig=figure();
    hold on;
    plot(Y(:,1),Y(:,3));
    plot(Y(:,7),Y(:,9));
    xlim([-1,2]);
    ylim([-1,1]);
    grid on;
    scatter([x0(1),x0(7)],[x0(3),x0(9)])
    %drawObsXZ(x0(7:9));
    hold off;
    saveas(fig, "D:\storage\CODE\matlab\framesXZ\orbitXZ_"+param.Wl+".png","png");

    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol,'Events', @inObsDE);
    tspan = [0, 1];
    [t1, Y1,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=TE;
    tspan = [0, -1];
    [t2, Y2,TE,YE,IE] = ode45(@Dfl, tspan, x0, optionsODE);
    T=T-TE;

    DifY=[flipud(Y2);Y1];
    Dift=[flipud(t2);t1];
    Dift=Dift-min(Dift);
    

    disp(["frames  :  ",length(DifY(:,1))])
    plot_inter_t(Dift,DifY,100/ut);
end

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
function T=plotSol(X,tmax)
    global param;
    %after :

    [Xl,Vl]=MoonPos(X(1));%toCart(param.al,param.el,param.wl,param.Wl,param.Il,X(1));
    %XP3=[param.DP3,0,0];
    %XP1=[param.DP1,0,0];
    Xts=SolPos(X(1));
    XP3=(param.DP3/param.Dter)*([Xl(1),Xl(2),Xl(3)]-Xts);
    XP1=(param.DP1/param.Dter)*([Xl(1),Xl(2),Xl(3)]-Xts);

    %radius:
    R=Xl+(1-X(2))*(XP3)+X(2)*XP1;
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*sin(X(4)),sin(X(4))*sin(X(3)),cos(X(4))]*Vn;

    x0=[R,V,Xl,Vl]';
    
    if(inObsDE(X(1),x0)<0)
        disp("probleme");
        T=0;
        return;
    end
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);
    tspan = [X(1), tmax];
    [t, Y] = ode45(@Dfl, tspan, x0, optionsODE);
    %disp("ooooooooooook");
    plot(Y(:,1),Y(:,2));

end

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
    Vn=X(5);%1*sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*sin(X(4)),sin(X(4))*sin(X(3)),cos(X(4))]*Vn;

    x0=[R,V]';
end

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

%==========================================================================
function [Opt,val]=optiFinder(fun,x0,dt,steps)
    Opt=x0;
    N=length(x0);
    val=fun(x0);
    for i= 0:steps
        %evaluate derivative
        der=zeros(0,N);
        for j=1:N
            Opt(j)=Opt(j)+1e-2;
            der(j)=(fun(Opt)-val)/1e-2;
            Opt(j)=Opt(j)-1e-2;
        end
        Opt=Opt+der*dt/norm(der);
        val2=fun(Opt);
        if(val2<val)
            Opt=Opt-der*dt/norm(der);
            dt=dt/2;
        else
            while(val2>val)
                val=val2;
                Opt=Opt+der*dt/norm(der);
                disp("value");
                disp(Opt);
                val2=fun(Opt);
                disp(der);
                disp(i);
                disp(dt);
                disp(val);
            end
            Opt=Opt-der*dt/norm(der);
        end
        disp("value");
        disp(Opt);
        disp(der);
        disp(i);
        disp(dt);
        disp(val);
        
    end
end

%==========================================================================
function [e,w,W,i,f]=getDeviceOrbital(X)
    global param;
    %after :

    Xl=toCart(param.al,param.el,param.wl,param.Wl,param.Il,X(1));
    %XP3=[param.DP3,0,0];
    %XP1=[param.DP1,0,0];
    
    XP3=(param.DP3/param.Dter)*[param.Dter+Xl(1),Xl(2),Xl(3)];
    XP1=(param.DP1/param.Dter)*[param.Dter+Xl(1),Xl(2),Xl(3)];

    %radius:
    R=Xl+(1-X(2))*(XP3)+X(2)*XP1;
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*cos(X(4)),cos(X(4))*sin(X(3)),sin(X(4))]*Vn;
    
    vr=dot(V,R)/Rn;
    vo=sqrt(Vn*Vn-vr*vr);

    h_vec=cross(R,V);
    h=norm(h_vec);


    i=acos(h_vec(3) / h);
    K = [0, 0, 1];
    N_vec = cross(K, h_vec);
    N = norm(N_vec);
    W = acos(N_vec(1)/N);
    if(N_vec(3)<0)
        W=-W;
    end

    e_vec = cross(V, h_vec) / param.Muter- R / Rn;
    e = norm(e_vec);

    w = acos(dot(N_vec, e_vec) / (N * e));
    if(e_vec(3)<0)
        w=-w;
    end

    f = acos(dot(e_vec, R) / (e * Rn));
    if(vr<0)
        f=-f;
    end
end
