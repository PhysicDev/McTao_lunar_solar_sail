
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
el=0.0;%54;
wl=0;%pi/2;
al=384399000;

%coeficient des fonctions :

% origine sur le point X3

p1=HP2/(DP2-DP3);
p2=HP2/(DP1-DP2);

off2=HP2+(DP2-DP3)*p2;

ud=al;
ut=2*pi*sqrt(al^3/Muter);

periodRatio=2;

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
param.nl=sqrt(param.Muter/param.al^3);
param.periodRatio=periodRatio;
param.as=param.al/(param.periodRatio^(2/3));
param.ns=sqrt(param.Muter/param.as^3);
param.timeStep=300/ut;
param.prec=1e-2/ut;


Xtest=[4*pi/6,0.5,205/360*2*pi];
%T=ObsTimeDE(Xtest)*ut
%T=ObsTime(Xtest)*ut


x0=[1.5*pi/3,0.5,206/360*2*pi];
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton',"Display","iter-detailed");
fun = @(x) -ObsTimeDE(x); % Negate the function to find maximum

[x, y] = fminunc(fun, x0, options);

disp(x);
disp(-y*ut);
disp(ObsTimeDE(x)*ut);

fig=figure();
plotOrbit(x);

saveas(fig, "D:\storage\CODE\matlab\orbit.png", "png");

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
    %after :

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

    x0=[R,V,Xl,getSpeed(param.al,param.el,param.wl,X(1))]';
    %options = odeset('Events', @inObsDE);
    reltol=1e-12;
    abstol=1e-12;
    
    optionsODE = odeset('RelTol', reltol, 'AbsTol', abstol);
    tspan = [0, 2];
    [t, Y] = ode45(@Df, tspan, x0, optionsODE);
    hold on;
    plot(Y(:,1),Y(:,2));
    grid on;

    %real orbit :
    f=0:0.1:(2*pi);

    R=param.as*(1-e*e)./(1+e*cos(f));
    plot(cos(f+w).*R,sin(f+w).*R);



    plot(Y(:,7),Y(:,8));


    R=param.al*(1-param.el*param.el)./(1+param.el*cos(f));
    plot(cos(f+param.wl).*R,sin(f+param.wl).*R);
    hold off;
end

function T=ObsTimeDE(X)
    global param;
    %after :

    Xl=toCart(param.al,param.el,param.wl,X(1));
    XP3=[param.DP3,0,0];
    XP1=[param.DP1,0,0];
    
    %radius:
    R=Xl+XP3*(1-X(2))+XP1*X(2);
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*Vn,sin(X(3))*Vn,0];

    x0=[R,V,Xl,getSpeed(param.al,param.el,param.wl,X(1))]';
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

%
% X : 
%
% position de la lune [0,2pi]
% position du satellite dans la zone [0,1]
% direction de la vitesse [0,2*pi]
%
function T=ObsTime(X)
    %disp(X) 8.9777    0.4037   -6.7428
    global param;
    %finding orbital component of orbit :
    %moonPos
    Xl=toCart(param.al,param.el,param.wl,X(1));
    XP3=[param.DP3,0,0];
    XP1=[param.DP1,0,0];
    Tl=perigeeTime(param.el,param.nl,X(1));
    
    %radius:
    R=Xl+XP3*(1-X(2))+XP1*X(2);
    Rn=norm(R);
    Vn=sqrt(param.Muter*(2/Rn-1/param.as));%norme de la vitesse fixé car demi grand axe fixé
    V=[cos(X(3))*Vn,sin(X(3))*Vn,0];
    

    h=cross(R,V);
    direction=sign(h(3));%used to know if we are moving against or in the same direction as the moon
    E=[(V(2)*h(3))/param.Muter-R(1)/Rn,(-V(1)*h(3))/param.Muter-R(2)/Rn,0];
    e=norm(E);
    w=atan2(E(2),E(1));
    f=acos(dot(E,R)/(e*Rn));
    Rangle=atan2(R(2),R(1));
    if(sign(sin(X(3)-Rangle)*cos(X(3)-Rangle))<0)
        f=-f;
    end
    
    Ts=perigeeTime(e,param.ns,f);

    Tup=0;
    Tdown=0;


    %initial Move (big step to find frontier with low precision) 
    %this step is used because its possible to have have multiple encounter zone , 
    %this way we avoid jumping between encouter zone in one step and having
    %a wrong observation time very long (which would be a maximum)
        %after:
        posS=R;
        posObs=Xl+XP3;
        while InObs(posS-posObs);
            Tup=Tup+param.timeStep;
            %updating positions
            posS=toCart(param.as,e,w,posInOrb(e,param.ns,direction*Tup+Ts));
            posObs=toCart(param.al,param.el,param.wl,posInOrb(param.el,param.nl,Tup+Tl))+XP3;
        end
        %before:
        posS=R;
        posObs=Xl+XP3;
        while InObs(posS-posObs);
            Tdown=Tdown-param.timeStep;
            %updating positions
            posS=toCart(param.as,e,w,posInOrb(e,param.ns,direction*Tdown+Ts));
            posObs=toCart(param.al,param.el,param.wl,posInOrb(param.el,param.nl,Tdown+Tl))+XP3;
        end

    %secondary Move (dichotomic search around frontier for high precision frontier)
        step=param.timeStep;
        while(step>param.prec)
            step=step/2;
            %after:
            posS=toCart(param.as,e,w,posInOrb(e,param.ns,direction*Tup+Ts));
            posObs=toCart(param.al,param.el,param.wl,posInOrb(param.el,param.nl,Tup+Tl))+XP3;
            %disp(Tup);
            %disp("pos")
            %disp(Ts+Tup-Ts);
            %disp(posInOrb(e,param.ns,direction*Tup+Ts));
            %disp(posS-posObs);
            %disp(R)
            if(InObs(posS-posObs))
                Tup=Tup+step;
            else
                Tup=Tup-step;
            end

            %before:
            posS=toCart(param.as,e,w,posInOrb(e,param.ns,direction*Tdown+Ts));
            posObs=toCart(param.al,param.el,param.wl,posInOrb(param.el,param.nl,Tdown+Tl))+XP3;
            if(InObs(posS-posObs))
                Tdown=Tdown-step;
            else
                Tdown=Tdown+step;
            end
        end
    %add up final time
    T=Tup-Tdown;

end