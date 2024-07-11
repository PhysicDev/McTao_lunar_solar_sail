
"""
value are normalized : 

distance unit : moon semi major axis
time unit : moon period (when considering the moon with a keplerian orbit)

"""



"""
loading library
"""

println("loading lib");
import SciMLBase
import OrdinaryDiffEq
println("opti")
using OptimalControl
using Plots
using CSV
using DataFrames

#
# gravity function : compute the graviational acceleration of an object
# X: pos object
# Pb: pos body
# mub: mu body
#
function gravity(X,Pb,mub)
    Rpos=Pb-X;
    Norm=sqrt(sum(Rpos.^2))
    return mub/(Norm^3)*Rpos;
end

#
# compute the derivative of an object subject to earth gravity
#
#
function lorenzVide!(du,u,p,t)
    du[1:3] = u[4:6];
    du[4:6] = gravity(u[1:3],[0,0,0],muTer);
end

#
# compute the derivative of an object subject to earth and moon gravity with an acceleration
#
# parameter p : (A,Mpos)
#
# A acceleration (function that return a vector of size 3)
# Mpos Moon data (function that return a vector of size 6)
#
function lorenz!(du,u,p,t)
    du[1:3] = u[4:6];
    du[4:6] = gravity(u[1:3],[0,0,0],muTer)+p[1](t) + gravity(u[1:3],p[2](t)[1:3],muLun)
end

#constant
#make this true if you want to solve the problem
resolve=true;
#make this true if you want to load the CSV file
load=true;
#make this true if you want to compute the moon dynamics
reloadMoon=true;

#loading CSV data
println("loading data");
if(load)
    times=CSV.File("./times.csv"; header=false) |> DataFrame #observations time info 
    startObservation=CSV.File("./startOBS.csv"; header=false) |> DataFrame #device position and speed when entering observation zone
    endObservation=CSV.File("./endOBS.csv"; header=false) |> DataFrame #device position and speed when exiting observation zone
    moonPos=CSV.File("./moonval.csv"; header=false) |> DataFrame #position and velocity of moon at time t=0
end

#initiate variable : 
startPoint=3 #starting observation ID
endPoint=5 #ending observation ID
Amax=5 #maximum acceleration (unused when minimizing Amax)
muTer=1.0*39.478417604357425 #acceleration constant of earth
muLun=1*0.483520433963301 #acceleration constant of moon

#compute the trajectory of the moon over 20 revolution (increase if needed)
println("compute moon trajectory")
if(reloadMoon)
    tspan = (0,20)
    prob3 = OrdinaryDiffEq.ODEProblem(lorenzVide!,moonPos[!,1],tspan)
    odeMoon3=OrdinaryDiffEq.solve(prob3,OrdinaryDiffEq.Tsit5(), reltol=1e-12, abstol=1e-12)
end

#
# on part de la sortie de la zone d'observation de départ
# jusqu'à l'entrée de la zone d'observation d'arrivée
#
ENDOBS=startObservation[!,endPoint]
STARTOBS=endObservation[!,startPoint]

tstart=times[!,startPoint][3]#0.353883781037970
tend=times[!,endPoint][2]#1.409239769931626


#objective value of the problem (0 mean not computed)
objective_output=zeros(24,24)


for j in range(1,22,step=1)
    for i in range(j+4,min(j+5,24),step=1)

        #start and end point of observation
        startPoint=j;
        endPoint=i;
        println(string("SP : ",startPoint," endP : ",endPoint));
        
        #gathering starting and end condition
        ENDOBS=startObservation[!,endPoint]
        STARTOBS=endObservation[!,startPoint]
        tstart=times[!,startPoint][3]
        tend=times[!,endPoint][2]
        
        #defining the problem
        @def ocp begin
            #tf ∈ R, variable
            maxia ∈ R, variable
            tf=tend-tstart
            ts=tstart-tstart
            #tend=0.728470473745049
            t ∈ [ ts, tf ], time
            x ∈ R^6, state
            u ∈ R^3, control
            #tf ≥ 0
            sum(u(t).^2)≤1
            10 ≥ maxia ≥ 0

            qx = x[1]
            qy = x[2]
            qz = x[3]
            vx = x[4]
            vy = x[5]
            vz = x[6]
        
            qx(ts) == STARTOBS[1]
            qy(ts) == STARTOBS[2]
            qz(ts) == STARTOBS[3]
        
            vx(ts) == STARTOBS[4]
            vy(ts) == STARTOBS[5]
            vz(ts) == STARTOBS[6]
        
            qx(tf) == ENDOBS[1]
            qy(tf) == ENDOBS[2]
            qz(tf) == ENDOBS[3]
        
            vx(tf) == ENDOBS[4]
            vy(tf) == ENDOBS[5]
            vz(tf) == ENDOBS[6]
        
            -10 ≤ qx(t) ≤ 10,      (1)
            -100 ≤ vx(t) ≤ 100,      (2)
            -10 ≤ qy(t) ≤ 10,      (3)
            -100 ≤ vy(t) ≤ 100,      (4)
            -10 ≤ qy(t) ≤ 10,      (5)
            -100 ≤ vy(t) ≤ 100,      (6)
            ẋ(t) == vcat([vx(t),vy(t),vz(t)],u(t)*maxia+gravity([qx(t),qy(t),qz(t)],[0,0,0],muTer)+ gravity(x(t)[1:3],odeMoon3(t+tstart)[1:3],muLun))
            #tf → min
            maxia → min
        end

        #sampling range
        sampling=[50,100,200,300];

        #solution and initial condition
        sol=nothing
        initg=nothing

        #ok 
        println("compute x0")
        tspan=[tstart,tend]
        u0(t)=[0,0,0];
        P=(u0,odeMoon3)

        #compute keplerian orbit
        tspan=[tstart,tend]
        pbtemp = OrdinaryDiffEq.ODEProblem(lorenzVide!,STARTOBS,tspan,P)
        sat1=OrdinaryDiffEq.solve(pbtemp,OrdinaryDiffEq.Tsit5(), reltol=1e-15, abstol=1e-15)
        
        tspan=[tend,tstart]
        pbtemp = OrdinaryDiffEq.ODEProblem(lorenzVide!,ENDOBS,tspan,P)
        sat2=OrdinaryDiffEq.solve(pbtemp,OrdinaryDiffEq.Tsit5(), reltol=1e-15, abstol=1e-15)

        # computing initial condition 
        # for now only computing dynamics without control 
        #
        probSat = OrdinaryDiffEq.ODEProblem(lorenz!,STARTOBS,tspan,P)
        x0=OrdinaryDiffEq.solve(probSat,OrdinaryDiffEq.Tsit5(), reltol=1e-15, abstol=1e-15)
        trueX0(t)=x0(t);

        trueX0(t)=(sat1(t)*(tend-t)+sat2(t)*(t-tstart))/(tend-tstart)
        initg = (state=trueX0, control=u0, variable=0)


        println("solve pb");#solving problem
        for samp in sampling 
            println(string("grid size : ",samp));
            if(resolve)
                if(isnothing(initg))#if there are no initial condition do nothing
                    sol = solve(ocp,display=false,grid_size=samp)
                else
                    sol = solve(ocp,display=false,grid_size=samp,init=initg)
                end
            end
            if(0==cmp(sol.message,"Solve_Succeeded"))#if we solve we set the new intitial solution then we continue
                println(string("solved  in ",sol.iterations," iterations || ",sol.objective))
                initg=(state=sol.state, control=sol.control, variable=sol.variable)
            else#if we don't solve we do nothing
                println("not solved :(")
            end
        end
        println(string("plot sol ... optimum : ",sol.objective));

        #plot(sol, size=(1200, 1800))
        objective_output[i,j]=sol.objective;#setting objective output variable
        #gui()



        #ploting solutions
        Tvec=[i for i in range(1,200)]/200*(tend-tstart);
        Tvec=[i for i in Tvec]

        states=[sol.state(x) for x in Tvec]

        X=[x[1] for x in states]
        Y=[x[2] for x in states]
        Z=[x[3] for x in states]

        mini=min(minimum(X),minimum(Y),minimum(Z))-0.1
        maxi=max(maximum(X),maximum(Y),maximum(Z))+0.1


        plot(X,Y,Z,label="trajectory",title=string("transfert from observation ",startPoint," to observation ",endPoint),xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi));
        scatter!([STARTOBS[1]],[STARTOBS[2]],[STARTOBS[3]],label="starting point")
        scatter!([ENDOBS[1]],[ENDOBS[2]],[ENDOBS[3]],label="final point")
        scatter!([0],[0],[0],label="Earth")

        #ploting start and end orbit:
        plot!(sat1,vars=(1,2,3),label="starting obrit")
        plot!(sat2,vars=(1,2,3),label="ending obrit")
        



        #title("transfert from observation "+startPoint+" to observation "+endPoint)
        savefig(string("images/Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));

        plot(sol, size=(1200, 1800))
        savefig(string("solGraph/Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));

        controlVal=[sol.control(t)*sol.objective for t in sol.times];
        X=[x[1] for x in controlVal]
        Y=[x[2] for x in controlVal]
        Z=[x[3] for x in controlVal]
        mini=-sol.objective*1.1;
        maxi=-mini;
        plot(X,Y,Z,label="control trajectory",xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi))
        scatter!([controlVal[1][1]],[controlVal[1][2]],[controlVal[1][3]],label="initial point")
        scatter!([controlVal[end][1]],[controlVal[end][2]],[controlVal[end][3]],label="end point")
        savefig(string("control/Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));

        
        controlVal=[sqrt(sum(sol.control(t).^2))*sol.objective for t in sol.times];
        plot(controlVal,ylim=(0,sol.objective*1.1));
        savefig(string("control/NORM_Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));

        
    end
end  
println("finished !!");
CSV.write("result.csv",  Tables.table(objective_output), writeheader=false)

gui()
while true
    sleep(0.1)
end
