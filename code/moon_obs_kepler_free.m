
Rsol=696000000;
Tsol=5777;
Rlun=1737000;
Rter=6500000;

Muter=6.7e-11*5.9737e24;

c=299792458;

%unité astronomique
UA=149597870690;

Dter=UA;

alpha=0.05;

DP1=Dter*Rlun/(Rsol-Rlun);

DP3=Dter*Rlun/(Rsol*(1+alpha)-Rlun);

A1=asin(Rlun/DP1);
A2=asin(Rlun/DP3);

DP2=(DP1*tan(A1)+DP3*tan(A2))/(tan(A1)+tan(A2));
HP2=tan(A1)*(DP1-DP2);

%orbite lunaire : 
el=0.0;%54;
wl=0;%pi/2;
al=384399000;

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
param.al=al/ud;
param.wl=wl;
param.el=el;
param.nl=sqrt(param.Muter/param.al^3);
param.periodRatio=periodRatio;
param.as=param.al/(param.periodRatio^(2/3));
param.ns=sqrt(param.Muter/param.as^3);
param.timeStep=300/ut;
param.prec=1e-2/ut;

x0=[120/360*pi*2,0.5,208/360*pi*2];
disp("ok");
disp(ObsTimeDE(x0));
disp("ok2");
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton',"Display","iter-detailed");
fun = @(x) -ObsTimeDE(x); % Negate the function to find maximum

[x, y] = fminunc(fun, x0, options);

disp(x);
disp(-y*ut);
disp(ObsTimeDE(x)*ut);
out=zeros(1,S-1);

plotOrbit(x);

%saveas(fig, "D:\storage\CODE\matlab\orbit.png", "png");

function P3=getXP3(Xl)
    global param;
    P3=(param.DP3/param.Dter)*[param.Dter+Xl(1),Xl(2),Xl(3)];%*sqrt(param.Dter*(param.Dter+2*Xl(1))))*[1,Xl(2)/(param.Dter+Xl(1)),0];%[1,X(2)/(X(1)+param.Dter),0];
end

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
        drawObs(Y(7:9)');
        grid on;

        title("T="+Ct*ut/3600+"h");
        hold off;
        saveas(fig, "D:\storage\CODE\matlab\framesXY\frame_"+i+".png","png");
        close(fig);

        Ct=Ct+step;
    end
end

function out=InObs(x)
    global param;
    d=sqrt(x(2)*x(2)+x(3)*x(3));
    out=all([d<param.p1*x(1) && -d<param.p1*x(1) && d<(param.off2-param.p2*x(1)) && -d<(param.off2-param.p2*x(1))]);
end

function [value, isterminal, direction]=inObsDE(t,X)
    global param;
    P3=(param.DP3/param.Dter)*[param.Dter+X(7),X(8),X(9)];%*sqrt(param.Dter*(param.Dter+2*Xl(1))))*[1,Xl(2)/(param.Dter+Xl(1)),0];%[1,X(2)/(X(1)+param.Dter),0];
    %P1=(param.DP1/param.Dter)*[param.Dter+X(7),Xl(8),Xl(9)];

    PX=dot(X(1:3)'-X(7:9)'-P3,P3)/norm(P3)^2*P3;
    PY=X(1:3)'-X(7:9)'-P3-PX;
    POS=[X(1)-X(7),X(2)-X(8),X(3)-X(9)]-P3;
    d=norm(PY);%sqrt(POS(2)*POS(2)+POS(3)*POS(3));%norm(POS());
    sx=norm(PX);%POS(1);
    out=all([d<param.p1*sx && -d<param.p1*sx && d<(param.off2-param.p2*sx) && -d<(param.off2-param.p2*sx)]);

    value=out-0.5;
    isterminal=1;
    direction = -1;
end

function V=getSpeed(a,e,w,f)
    global param;
    cst=sqrt(param.Muter/(a*(1-e*e)));
    V=cst*[-sin(w+f)-e*sin(w),cos(w+f)+e*cos(w),0];
end

function Y=Df(t,X)
    global param;
    accS=param.Muter/norm(X(1:3))^3;
    accL=param.Muter/norm(X(7:9))^3;
    Y=[X(4:6);-accS*X(1:3);X(10:12);-accL*X(7:9)];
end


function out=plotOrbit(X)
    global param;
    global ut;
    %after :

    Xl=toCart(param.al,param.el,param.wl,X(1));
    %XP3=[param.DP3,0,0];
    %XP1=[param.DP1,0,0];
    
    XP3=(param.DP3/param.Dter)*[param.Dter+Xl(1),Xl(2),Xl(3)];
    XP1=(param.DP1/param.Dter)*[param.Dter+Xl(1),Xl(2),Xl(3)];

    %radius:
    R=Xl+(1-X(2))*(XP3)+X(2)*XP1;
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3)),sin(X(3)),0]*Vn;

    x0=[R,V,Xl,getSpeed(param.al,param.el,param.wl,X(1))]';
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    %fig=figure();
    fig=figure('Visible', 'off');
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);
    tspan = [0, 2];
    [t, Y] = ode45(@Df, tspan, x0, optionsODE);
    hold on;
    plot(Y(:,1),Y(:,2));
    plot(Y(:,7),Y(:,8));
    xlim([-1,2]);
    ylim([-1,1]);
    grid on;
    scatter([x0(1),x0(7)],[x0(2),x0(8)])
    %drawObs(x0(7:9));
    hold off;
    %saveas(fig, "D:\storage\CODE\matlab\framesXY\orbitXY_"+param.Wl+".png","png");
    

    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol,'Events', @inObsDE);
    tspan = [0, 1];
    [t1, Y1,TE,YE,IE] = ode45(@Df, tspan, x0, optionsODE);
    T=TE;
    tspan = [0, -1];
    [t2, Y2,TE,YE,IE] = ode45(@Df, tspan, x0, optionsODE);
    T=T-TE;

    DifY=[flipud(Y2);Y1];
    Dift=[flipud(t2);t1];
    Dift=Dift-min(Dift);
    

    %disp(["frames  :  ",length(DifY(:,1))])
    
    
    plot_inter_t(Dift,DifY,100/ut);
    
end

function drawObs(Xl)
    %on desine en 2D donc on neglige l'axe z pour l'instant
    global param;
    Xls=[param.Dter+Xl(1),Xl(2),Xl(3)];
    XP3=(param.DP3/param.Dter)*Xls;
    XP1=(param.DP1/param.Dter)*Xls;
    XP2x=(param.DP2/param.Dter)*Xls;
    XP2y=(XP2x/norm(XP2x))*[0,-1,0;1,0,0;0,0,1]*param.HP2;
    poly=[XP1+Xl';Xl'+XP2x+XP2y;Xl'+XP3;Xl'+XP2x-XP2y;XP1+Xl'];
    plot(poly(:,1),poly(:,2));


    xlim([min(poly(:,1)),max(poly(:,1))])
    ylim([min(poly(:,2)),max(poly(:,2))])
end

function T=ObsTimeDE(X)
    global param;
    %after :

    Xl=toCart(param.al,param.el,param.wl,X(1));
    %XP3=[param.DP3,0,0];
    %XP1=[param.DP1,0,0];
    Vl=getSpeed(param.al,param.el,param.wl,X(1));
    
    XP3=(param.DP3/param.Dter)*[param.Dter+Xl(1),Xl(2),Xl(3)];
    XP1=(param.DP1/param.Dter)*[param.Dter+Xl(1),Xl(2),Xl(3)];

    %radius:
    R=Xl+(1-X(2))*(XP3)+X(2)*XP1;
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*Vn,sin(X(3))*Vn,0];

    x0=[R,V,Xl,Vl]';
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol,'Events', @inObsDE);
    tspan = [0, 1];
    [t, Y,TE,YE,IE] = ode45(@Df, tspan, x0, optionsODE);
    T=TE;
    tspan = [0, -1];
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

function X=toCart(a,e,w,f)
    if(e==1)
        e=0.9999999;%pour eviter les div par zéro
    end
    r=a*(1-e*e)/(1+e*cos(f));
    x=cos(f+w)*r;
    y=sin(f+w)*r;
    z=0;%on neglige la 3eme dimension
    X=[x,y,z];
end

function [e,w,f]=getDeviceOrbital(X)
    global param;

    Xl=toCart(param.al,param.el,param.wl,X(1));
    XP3=[param.DP3,0,0];
    XP1=[param.DP1,0,0];

    %radius:
    R=Xl+XP3*(1-X(2))+XP1*X(2);
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*Vn,sin(X(3))*Vn,0];
    

    h=cross(R,V);
    E=[(V(2)*h(3))/param.Muter-R(1)/Rn,(-V(1)*h(3))/param.Muter-R(2)/Rn,0];
    e=norm(E);
    w=atan2(E(2),E(1));
    f=acos(dot(E,R)/(e*Rn));

end
