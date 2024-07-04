println("loading lib");

import SciMLBase
import OrdinaryDiffEq
using OptimalControl
using Plots
using CSV
using DataFrames

#
# X: pos object
# Pb: pos body
# mub: mu body
#
function gravity(X,Pb,mub)
    Rpos=Pb-X;
    Norm=sqrt(sum(Rpos.^2))
    return mub/(Norm^3)*Rpos;
end
function lorenzVide!(du,u,p,t)
    du[1:3] = u[4:6];
    du[4:6] = gravity(u[1:3],[0,0,0],muTer);
end


resolve=true;
load=true;
reloadMoon=true;
println("loading data");
if(load)
    times=CSV.File("./times.csv"; header=false) |> DataFrame
    startObservation=CSV.File("./startOBS.csv"; header=false) |> DataFrame
    endObservation=CSV.File("./endOBS.csv"; header=false) |> DataFrame
    moonPos=CSV.File("./moonval.csv"; header=false) |> DataFrame
end

startPoint=3
endPoint=5
Amax=5
muTer=1.0*39.478417604357425
muLun=1*0.483520433963301

println("compute moon trajectory")
if(reloadMoon)
    tspan = (0,20)
    prob3 = OrdinaryDiffEq.ODEProblem(lorenzVide!,moonPos[!,1],tspan)
    odeMoon3=OrdinaryDiffEq.solve(prob3,OrdinaryDiffEq.Tsit5(), reltol=1e-12, abstol=1e-12)
end
println("define pb");

#
# on part de la sortie de la zone d'observation de départ
# jusqu'à l'entrée de la zone d'observation d'arrivée
#
ENDOBS=startObservation[!,endPoint]


"""[-0.167269293387500,
        1.09832700578628,
        0.0450395299341377,
        -3.13998323037230,
        -4.62822855823360,
        -0.450462273954064];"""

STARTOBS=endObservation[!,startPoint]


"""[0.231357458212836,
          0.922293662828187,
          0.0670351958926985,
          -5.50343680893103,
          -3.61488252858483,
          -0.343544167374632]"""

tstart=times[!,startPoint][3]#0.353883781037970
tend=times[!,endPoint][2]#1.409239769931626


function F0(X)
    return vcat(X(4:6),gravity(X(1:3),[0,0,0],muTer))
end




for j in range(1,22,step=1)
    for i in range(j+2,min(j+5,24),step=1)
        startPoint=j;
        endPoint=i;
        println(string("SP : ",startPoint," endP : ",endPoint));
        
        ENDOBS=startObservation[!,endPoint]
        STARTOBS=endObservation[!,startPoint]
        tstart=times[!,startPoint][3]
        tend=times[!,endPoint][2]

        @def ocp begin
            #tf ∈ R, variable
            tf=tend
            #tend=0.728470473745049
            
            t ∈ [ tstart, tf ], time
            x ∈ R^6, state
            u ∈ R^3, control
            #tf ≥ 0
            sum(u(t).^2) ≤ Amax^2
            qx = x[1]
            qy = x[2]
            qz = x[3]
            vx = x[4]
            vy = x[5]
            vz = x[6]
        
            qx(tstart) == STARTOBS[1]
            qy(tstart) == STARTOBS[2]
            qz(tstart) == STARTOBS[3]
        
            vx(tstart) == STARTOBS[4]
            vy(tstart) == STARTOBS[5]
            vz(tstart) == STARTOBS[6]
        
            qx(tf) == ENDOBS[1]
            qy(tf) == ENDOBS[2]
            qz(tf) == ENDOBS[3]
        
            vx(tf) == ENDOBS[4]
            vy(tf) == ENDOBS[5]
            vz(tf) == ENDOBS[6]
        
            -999 ≤ qx(t) ≤ 999,      (1)
            -999 ≤ vx(t) ≤ 999,      (2)
            -999 ≤ qy(t) ≤ 999,      (3)
            -999 ≤ vy(t) ≤ 999,      (4)
            -999 ≤ qy(t) ≤ 999,      (5)
            -999 ≤ vy(t) ≤ 999,      (6)
            ẋ(t) == vcat([vx(t),vy(t),vz(t)],u(t)+gravity([qx(t),qy(t),qz(t)],[0,0,0],muTer)+ gravity(x(t)[1:3],odeMoon3(t)[1:3],muLun))
            #tf → min
            ∫( 0.5sum(u(t).^2) ) → min
        end

        println("solve pb");
        sampling=[200,300,400,500,600];
        sol=nothing
        for samp in sampling
            println(string("atempting grid size : ",samp));
            if(resolve)
                sol = solve(ocp,display=false,grid_size=samp)
            end
            if(0==cmp(sol.message,"Solve_Succeeded"))
                println(string("solved  in ",sol.iterations," iterations"))
                break
            end
            println("not solved :(")
        end
        println("plot sol");

        #plot(sol, size=(1200, 1800))

        #gui()

        Tvec=[i for i in range(1,200)]/200*(tend-tstart);
        Tvec=[i+tstart for i in Tvec]

        states=[sol.state(x) for x in Tvec]

        X=[x[1] for x in states]
        Y=[x[2] for x in states]
        Z=[x[3] for x in states]

        plot(X,Y,Z,label="trajectory",title=string("transfert from observation ",startPoint," to observation ",endPoint));
        scatter!([STARTOBS[1]],[STARTOBS[2]],[STARTOBS[3]],label="starting point")
        scatter!([ENDOBS[1]],[ENDOBS[2]],[ENDOBS[3]],label="final point")
        scatter!([0],[0],[0],label="Earth")
        #title("transfert from observation "+startPoint+" to observation "+endPoint)
        savefig(string("images/Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));

        plot(sol, size=(1200, 1800))
        savefig(string("solGraph/Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));
    end
end  
println("finished !!");


gui()
while true
    sleep(0.1)
end
