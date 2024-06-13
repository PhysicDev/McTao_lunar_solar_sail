
Rsol=696000000;
Tsol=5777;
Rlun=1737000;
Rter=6378000;

Muter=6.7e-11*5.9737e24;

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
param.periodSol=periodSol/ut;
param.Rter=Rter/ud;
param.Rlun=Rlun/ud;




%on opti !!!
Synperiod=1.082;
nbOpt=14;

Opt=zeros(10,nbOpt);

options = optimoptions('fminunc', 'Algorithm', 'quasi-newton',"Display","off","MaxFunctionEvaluations",1e5);
fun = @(x) BestVal(x); % Negate the function to find maximum
x0=0.1;
[t1, y1] = fminunc(fun, x0, options);
pos=0;
for i=1:(nbOpt/2)
    %opt1
    x0=t1+pos;
    [tn, yn] = fminunc(fun, x0, options);
    [out,Pm,Pz,V]=BestVal(tn);
    Opt(1,2*i-1)=tn;
    Opt(2:4,2*i-1)=Pm;
    Opt(5:7,2*i-1)=Pz;
    Opt(8:10,2*i-1)=V;

    pos=pos+Synperiod;

    %opt2
    x0=pos-t1;
    [tn, yn] = fminunc(fun, x0, options);
    [out,Pm,Pz,V]=BestVal(tn);
    Opt(1,2*i)=tn;
    Opt(2:4,2*i)=Pm;
    Opt(5:7,2*i)=Pz;
    Opt(8:10,2*i)=V;
end

DV=zeros(nbOpt,nbOpt);

Opt

DV;
I=1
J=3

for I=1:(nbOpt-1)
    for J=(I+1):min(nbOpt,I+5)
        I
        J
        %creating mesh
        S1 = Opt(1,I);
        E1 = Opt(1,J);
        S2 = Opt(1,I);
        E2 = Opt(1,J);
        
        prec=2e-2;
        
        Sdat = S1:prec:E1;
        Edat = S2:prec:E2;
        
        [Sgrid, Egrid] = meshgrid(Sdat, Edat);
        
        
        [Vi,Vf]=solvelambert(Opt(5:7,I)',Opt(5:7,J)',Opt(1,J)-Opt(1,I),-1,param.Muter);
        t=Opt(1,J)-Opt(1,I);
        Dv=norm(Vi-Opt(8:10,I))+norm(Vf-Opt(8:10,J));
        DV(I,J)=Dv;
       
        %plotOrbitXY(a1,e1,w1,W1,i1);
        optionsODE = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
        
        
        x0=[Opt(5:7,I)',Opt(8:10,I)'];
        tspan=[0,Opt(1,J)-Opt(1,I)];
        [time1, Y1] = ode45(@Df, tspan, x0,optionsODE);
        
        sol1 = ode45(@Df, tspan, x0,optionsODE);%struct('x', time1, 'y', Y1.');
        
        
        
        
        x0=[Opt(5:7,J)',Opt(8:10,J)'];
        tspan=[0,Opt(1,I)-Opt(1,J)];
        [time2, Y2] = ode45(@Df, tspan, x0,optionsODE);
        sol2 = ode45(@Df, tspan, x0,optionsODE);%struct('x', time2, 'y', Y2.');
        
        num_departure_dates = length(Sdat);
        num_arrival_dates = length(Edat);
        delta_v = zeros(num_departure_dates, num_arrival_dates);
        
        dt=[];
        Ps1=[];
        Ps2=[];
        Vs2=[];
        Vs3=[];
        Px=[];
        Py=[];
        for pci = 1:(num_departure_dates-3)
            departure_date = Sdat(pci);
            int=pci+3:num_departure_dates;

            Py=[Py,int];
            Px=[Px;repmat(pci,length(int),1)];
        
            subEdat=Edat(int);
            Ys2=deval(sol2,subEdat-Opt(1,J))';
            Ys1=repmat(deval(sol1,departure_date-Opt(1,I))',length(subEdat),1);
            %Ys1=deval(sol1,departure_date-Opt(1,I))'
        
            Ps1=[Ps1;Ys1(:,1:3)];
            Vs1=[Vs1;Ys1(:,4:6)];
            
            Ps2=[Ps2;Ys2(:,1:3)];
            Vs2=[Vs2;Ys2(:,4:6)];


            dt=[dt;(subEdat-departure_date)'];
        end

        

        [Vi,Vf]=solvelambert(Ps1,Ps2,dt,ones(length(dt), 1),param.Muter);%= compute_delta_v(departure_date, arrival_date, departure_body, arrival_body, mu_sun);
        
        

        for pos = 1:length(Px)
            delta_v(Px(pos),Py(pos))=vecnorm((Vi(pos,:)-Vs1(pos,:))')+vecnorm((Vf(pos,:)-Vs2(pos,:))');
        end
        
        delta_v=delta_v*ud/ut/1000;
        
        f=figure("Visible","off");
        contourf(Sgrid,Egrid, delta_v,(1:50)/10, 'LineColor', 'none');
        colorbar;
        xlabel('Arrival Date');
        ylabel('Departure Date');
        title('Porkchop Plot');
        saveas(f, "D:\storage\CODE\github\McTao_lunar_solar_sail\images\porkchop\porkchop2_"+I+"_"+J+".png", "png");
        close(f);
        delta_v(delta_v==0)=1000;
        DV(I,J)=min(delta_v(:));
        %plot(Y(:,1),Y(:,2));
        %hold on;
        %scatter([Opt(5,I),Opt(5,J)],[Opt(6,I),Opt(6,J)])
        %hold off;
    end
end

function T=perigeeTime(e,n,f)
    E=2*atan(sqrt((1-e)/(1+e))*tan(f/2));
    M=E-e*sin(E);
    T=M/n;
end

function f=posInOrb(e,n,t)
    M=t*n;
    options = optimoptions('fsolve', ...
                            'Display','off', ...
                            'FunctionTolerance',1e-8);
    E=fsolve(@(x) (x-e*sin(x)-M),M,options);
    f=2*atan(sqrt((1+e)/(1-e))*tan(E/2));%sqrt((1+e)/(1-e))
end


function Y=Df(t,X)
    global param;
    accL=param.Muter/norm(X(1:3))^3;
    Y=[X(4:6);-accL*X(1:3)];
end

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

    x=r.*(CW*CT-SW*ST*CI);
    y=r.*(SW*CT+CW*ST*CI);
    z=r.*ST*sin(i);%on neglige la 3eme dimension
    X=[x',y',z'];
end

function drawObsXZ(Xl)
    %on desine en 2D donc on neglige l'axe z pour l'instant
    global param;
    Xls=[param.Dter+Xl(1),Xl(2),Xl(3)];
    XP3=(param.DP3/param.Dter)*Xls;
    XP1=(param.DP1/param.Dter)*Xls;
    XP2x=(param.DP2/param.Dter)*Xls;
    XP2y=(XP2x/norm(XP2x))*[0,0,-1;0,1,0;1,0,0]*param.HP2;
    poly=[XP1+Xl';Xl'+XP2x+XP2y;Xl'+XP3;Xl'+XP2x-XP2y;XP1+Xl'];
    plot(poly(:,1),poly(:,3));


end

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

end

function [out,Pm,Pz,V]=BestVal(t)
    global param;
    
    %pos du soleil
    A=t/param.periodSol*2*pi+pi;
    Xst=[cos(A),sin(A),0]*param.Dter;

    if(t<0)
        t=1-mod(t,1);
    else
        t=mod(t,1);
    end
    %step 1 : calculer la position de la Lune (dur)
    T=perigeeTime(param.el,param.nl,param.fl);
    f=posInOrb(param.el,param.nl,t-T);
    posLun=toCart(param.al,param.el,param.wl,param.Wl,param.Il,f);

    %step 2: calculer la position de la zone d'observation (facile)
    P2=posLun+param.DP2/param.Dter*(Xst+posLun);
    %graphe:
    %plotMoon(t);
    %drawObs(posLun',Xst);
    %xlim([-2,2])
    %ylim([-2,2])
    %hold off;

    %retourner le payoff (très facile)
    out=(norm(posLun)-norm(P2))^2;
    Pm=posLun;
    Pz=P2;
    V=getSpeed(param.al,param.el,param.wl,param.Wl,param.Il,f);
end

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

function plotOrbitXY(a,e,w,W,i)
    angle=-pi:0.01:pi;
    X=toCart(a,e,w,W,i,angle);
    plot(X(:,1),X(:,2));
end

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
