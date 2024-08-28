
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
using NLPModelsIpopt
using Base.Filesystem
using JLD2, JSON3

function integral(f,t0,tf,N)
    h=(tf-t0)/N
    Tvec=[i-0.5 for i in range(1,N)]*h;
    Tvec=[t+t0 for t in Tvec]

    val=[f(t)*h for t in Tvec]
    return sum(val)
end

function AddFolder(folder_path::String)
    if !isdir(folder_path)
        try
            mkdir(folder_path)
            println("Directory created: $folder_path")
        catch e
            println("Failed to create directory: $folder_path. Error: $e")
        end
    else
        println("Directory already exists: $folder_path")
    end
end

AddFolder("images")
AddFolder("control")
AddFolder("solGraph")
AddFolder("solutions")

function integral(f,t0,tf,N)
    h=(tf-t0)/N
    Tvec=[i-0.5 for i in range(1,N)]*h;
    Tvec=[t+t0 for t in Tvec]

    val=[f(t)*h for t in Tvec]
    return sum(val)
end

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


function lorenzSun!(du,u,p,t)
    du[1:3] = u[4:6];
    du[4:6] = gravity(u[1:3],[0,0,0],muSol);
    du[4:6] = gravity(u[1:3],u[7:9],muLun);

    du[7:9] = u[10:12];
    du[10:12] = gravity(u[7:9],[0,0,0],muSol);
    du[10:12] = gravity(u[7:9],u[1:3],muTer);
end

function Mpos(sunSys,t)
    val=sunSys(t)
    orig=(val[1:6]*muTer+val[7:12]*muLun)/(muTer+muLun)
    return val[7:12]-orig
end

function Epos(sunSys,t)
    val=sunSys(t)
    orig=(val[1:6]*muTer+val[7:12]*muLun)/(muTer+muLun)
    return val[1:6]-orig
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
    val=p[2](t)
    orig=(val[1:6]*muTer+val[7:12]*muLun)/(muTer+muLun)
    Epos=val[1:6]-orig
    Lpos=val[7:12]-orig
    du[4:6] = gravity(u[1:3],Epos[1:3],muTer)+p[1](t) + gravity(u[1:3],Lpos[1:3],muLun)
end

function massDer(u)
    Norm=sqrt(sum(u.^2))
    return Norm
end


function acceleration(t,Pos,sunSys)
    val=sunSys(t)
    orig=(val[1:6]*muTer+val[7:12]*muLun)/(muTer+muLun)
    Epos=val[1:6]-orig
    Lpos=val[7:12]-orig
    return gravity(Pos,Epos[1:3],muTer)+ gravity(Pos,Lpos[1:3],muLun)
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
    moonEarthPos=CSV.File("./Sun_Moon_Earth_jan_2025.csv"; header=true) |> DataFrame #position and velocity of moon and Earth at time t=0 (2025)
    N=ncol(times);
end

#initiate variable : 
muTer=1*39.478417604357425 #acceleration constant of earth
muLun=1*0.483520433963301 #acceleration constant of moon
muSol=1*1.314418294336131e+07 #acceleration constant of sun

#device related variable (use data of RIT 10 evo thruster):
Isp=3000 #isp in s
Mass=50 #payload mass in kg
FuelMass=40 #fuel mass in kg
Thrust=0.025 #thrust force in N

#normalisation constant
ud=384399000 #distance : Earth Moon mean distance
ut=2.371834349258416e+06 #time Lunar mean period

#thruster related computed variable
g0=9.80665; 
Q=Thrust/Isp/g0;
NormThrust=0.02/ud*ut*ut;
NormQ=Q*ut;
println("normalized Thrust")
println(NormThrust)
println("flow : ")
println(Q)
println("normalized flow : ")
println(NormQ)

#compute the trajectory of the moon over 20 revolution (increase if needed)
println("compute moon trajectory")

if(reloadMoon)
    tspan = (0,times[!,end][3]+1)
    prob3 = OrdinaryDiffEq.ODEProblem(lorenzSun!,vcat(Array(moonEarthPos[3,:]),Array(moonEarthPos[2,:])),tspan)
    odeMoonEarth3=OrdinaryDiffEq.solve(prob3,OrdinaryDiffEq.Tsit5(), reltol=1e-12, abstol=1e-12)
end

#
# on part de la sortie de la zone d'observation de départ
# jusqu'à l'entrée de la zone d'observation d'arrivée
#

startPoint=3 #starting observation ID
endPoint=5 #ending observation ID

ENDOBS=startObservation[!,endPoint]
STARTOBS=endObservation[!,startPoint]

tstart=times[!,startPoint][3]#0.353883781037970
tend=times[!,endPoint][2]#1.409239769931626

fun(t)=t

if(resolve)
    #objective value of the problem (0 mean not computed)
    objective_output=zeros(N,N)
    deltas_V=zeros(N,N)
    mass=zeros(2,N)
    
    sols=Array{Any}(undef, N,N)
    #s = readline()
    for j in range(1,N-2,step=1)
        for k in [2]

            global ENDOBS,STARTOBS,startPoint,endPoint,tstart,fun,tend,Tvec,mini,maxi,tspan

            #start and end point of observation
            startPoint=j;
            i=j+k;
            endPoint=i;
            Mi=FuelMass+Mass;
            
            Pi=j%2+1;

            if(j>2)
                Mi=mass[Pi,j-2];
            end

            println("    ")
            println(string("startP : ",startPoint," endP : ",endPoint));
            
            #gathering starting and end condition
            ENDOBS=startObservation[!,endPoint]
            STARTOBS=endObservation[!,startPoint]
            tstart=times[!,startPoint][3]
            tend=times[!,endPoint][2]
            
            time0=tstart#-tstart
            timef=tend#-tstart
            #defining the problem
            @def ocp begin
                #tf ∈ R, variable
                #maxia ∈ R, variable

                tf=timef
                ts=time0

                #tend=0.728470473745049
                t ∈ [ ts, tf ], time
                x ∈ R^7, state
                u ∈ R^3, control
                #tf ≥ 0
                sum(u(t).^2)≤1
                #1 ≥ maxia ≥ 0

                rx = x[1]
                ry = x[2]
                rz = x[3]

                vx = x[4]
                vy = x[5]
                vz = x[6]

                m = x[7]
            
                rx(ts) == STARTOBS[1]
                ry(ts) == STARTOBS[2]
                rz(ts) == STARTOBS[3]
                
                vx(ts) == STARTOBS[4]
                vy(ts) == STARTOBS[5]
                vz(ts) == STARTOBS[6]

                m(ts) == Mi;
            
                rx(tf) == ENDOBS[1]
                ry(tf) == ENDOBS[2]
                rz(tf) == ENDOBS[3]
            
                vx(tf) == ENDOBS[4]
                vy(tf) == ENDOBS[5]
                vz(tf) == ENDOBS[6]
            
                -10 ≤ rx(t) ≤ 10,      (1)
                -100 ≤ vx(t) ≤ 100,      (2)
                -10 ≤ ry(t) ≤ 10,      (3)
                -100 ≤ vy(t) ≤ 100,      (4)
                -10 ≤ ry(t) ≤ 10,      (5)
                -100 ≤ vy(t) ≤ 100,      (6)
                Mass ≤ m(t) ≤ Mi,      (7)


                ẋ(t) == vcat( [vx(t),vy(t),vz(t)] #position derivative
                              , u(t)*NormThrust/m(t) #control acceleration
                              + acceleration(t+tstart-time0,[rx(t),ry(t),rz(t)],odeMoonEarth3) #gravity of Earth and Moon (see accelerationmethod for more info)
                              , [-NormQ*(FuelMass+Mass)/NormThrust * (max(0,(sum(u(t).^2))))^(1/2) ]
                            )
                #tf → min
                ∫( 0.5*(sum(u(t).^2)) ) → min
            end


            
            #sampling range
            sampling=[50,100,200];

            #solution and initial condition
            sol=nothing
            initg=nothing

            #ok 
            println("compute x0")
            tspan=[tstart,tend]
            u0(t)=[0,0,0];
            P=(u0,odeMoonEarth3)

            #compute keplerian orbit
            tspan=[tstart,tend]
            pbtemp = OrdinaryDiffEq.ODEProblem(lorenz!,STARTOBS,tspan,P)
            sat1=OrdinaryDiffEq.solve(pbtemp,OrdinaryDiffEq.Tsit5(), reltol=1e-15, abstol=1e-15)
            
            tspan=[tend,tstart]
            pbtemp = OrdinaryDiffEq.ODEProblem(lorenz!,ENDOBS,tspan,P)
            sat2=OrdinaryDiffEq.solve(pbtemp,OrdinaryDiffEq.Tsit5(), reltol=1e-15, abstol=1e-15)

            # computing initial condition 
            # for now only computing dynamics without control 
            #
            probSat = OrdinaryDiffEq.ODEProblem(lorenz!,STARTOBS,tspan,P)
            x0=OrdinaryDiffEq.solve(probSat,OrdinaryDiffEq.Tsit5(), reltol=1e-15, abstol=1e-15)
            trueX0(t)=x0(t);

            trueX0(t)=vcat((sat1(t+tstart-time0)*(tend-t+tstart-time0)+sat2(t+tstart-time0)*(t-time0))/(tend-tstart),[Mi])
            initg = (state=trueX0, control=u0)
            #initg=nothing
            println("    ");
            println("solve time pb");#solving problem

            #solving process
            for samp in sampling 
                println(string("grid size : ",samp));
                if(isnothing(initg))#if there are no initial condition do nothing
                    sol = solve(ocp,display=false,grid_size=samp,max_iter=3000)
                else
                    sol = solve(ocp,display=false,grid_size=samp,init=initg)
                end
                if(0==cmp(sol.message,"Solve_Succeeded"))#if we solve we set the new intitial solution then we continue
                    println(string("solved in ",sol.iterations," iterations || fuel used : ",Mi-sol.state(tend)[7] , "Kg : ",(Mi-sol.state(tend)[7])/FuelMass*100," %"))
                    initg=(state=sol.state, control=sol.control, variable=sol.variable)#we update the initial condition for the next solving iteration
                else#if we don't solve we do nothing
                    println("not solved :(")
                end
            end

            #printint result
            println(" ");
            println(string("founded solution max acceleration ",sol.objective));

            
            fun(t)=sqrt(sum((sol.control(t)).^2))/sol.state(t)[7];
            tempDV=integral(fun,time0,timef,10000)*NormThrust;
            println(string("delta V :",tempDV))
            deltas_V[i,j]=tempDV;
            println(string("plot sol ... optimum : ",sol.objective));

            mass[Pi,j]=sol.state(tend)[7];

            objective_output[i,j]=sol.objective;
           
            sols[j,i]=sol;

            save(sol, filename_prefix=string("./solutions/sol_",startPoint,"_",endPoint ))



            #ploting solutions
            Tvec=[i for i in range(1,1000)]/1000*(tend-tstart);
            LowTvec=[i for i in range(1,100)]/100*(tend-tstart);
            Tvec=[i+time0 for i in Tvec]
            LowTvec=[i+time0 for i in LowTvec]

            states=[sol.state(x) for x in Tvec]
            sols[j,i]=states;

            X=[x[1] for x in states]
            Y=[x[2] for x in states]
            Z=[x[3] for x in states]

            mini=min(minimum(X),minimum(Y),minimum(Z))*1.5
            maxi=max(maximum(X),maximum(Y),maximum(Z))*1.5

            massDat=[x[7] for x in states]
            plot(Tvec,massDat,title="fuel tank");
            savefig(string("images/mass_progress_",startPoint,"_",endPoint,".png"));
            #gui();
            #wait = readline()

            plot(X,Y,Z,label="trajectory",title=string("transfert from observation ",startPoint," to observation ",endPoint),xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi),size=(1920/2,1080/2));
            
            for t in LowTvec
                St=sol.state(t);
                Ct=sol.control(t);
                plot!([St[1],St[1]+Ct[1]],[St[2],St[2]+Ct[2]],[St[3],St[3]+Ct[3]],label=false,color=:red)
            end
            scatter!([STARTOBS[1]],[STARTOBS[2]],[STARTOBS[3]],label="starting point")
            scatter!([ENDOBS[1]],[ENDOBS[2]],[ENDOBS[3]],label="final point")
            scatter!([0],[0],[0],label="Earth")


            #ploting start and end orbit:
            plot!(sat1,vars=(1,2,3),label="starting obrit",color=:green)
            plot!(sat2,vars=(1,2,3),label="ending obrit",color=:orange)
            


            #title("transfert from observation "+startPoint+" to observation "+endPoint)
            savefig(string("images/Transfert_",startPoint,"_",endPoint,".png"));

            plot(sol, size=(1200, 1800))
            savefig(string("solGraph/Transfert_",startPoint,"_",endPoint,".png"));

            controlVal=[sol.control(t)*sol.objective for t in Tvec];
            
            Xc=[x[1] for x in controlVal]
            Yc=[x[2] for x in controlVal]
            Zc=[x[3] for x in controlVal]

            newX=[(states[i][4]*Xc[i]+states[i][5]*Yc[i]+states[i][6]*Zc[i])/sqrt(sum(states[i][4:6].^2)) for i in range(1,length(states))];
            Xax=[x[2]*x[6]-x[3]*x[5] for x in states]
            Yax=[x[3]*x[4]-x[1]*x[6] for x in states]
            Zax=[x[1]*x[5]-x[2]*x[4] for x in states]
            newZ=[(Xax[i]*Xc[i]+Yax[i]*Yc[i]+Zax[6]*Zc[i])/sqrt(Xax[i]^2+Yax[i]^2+Zax[i]^2) for i in range(1,length(states))];

            Xaxy=[Yax[i]*states[i][6]-Zax[i]*states[i][5]  for i in range(1,length(states))]
            Yaxy=[Zax[i]*states[i][4]-Xax[i]*states[i][6]  for i in range(1,length(states))]
            Zaxy=[Xax[i]*states[i][5]-Yax[i]*states[i][2]  for i in range(1,length(states))]
            
            newY=[(Xaxy[i]*Xc[i]+Yaxy[i]*Yc[i]+Zaxy[6]*Zc[i])/sqrt(Xaxy[i]^2+Yaxy[i]^2+Zaxy[i]^2) for i in range(1,length(states))];
            #newY=[sqrt(sum( ([Xc[i],Yc[i],Zc[i]] - newZ[i]*[Xax[i],Xax[i],Yax[i]]/sqrt(Xax[i]^2+Yax[i]^2+Zax[i]^2) - newX[i]*states[i][4:6]/sqrt(sum(states[i][4:6].^2))).^2 )) for i in range(1,length(states))];
            mini=-sol.objective*1.1;
            maxi=-mini;
            #println(newY)
            plot(newX,newY,newZ,label="control trajectory",xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi))
            #println([newX[1]],[newY[1]],[newZ[1]])
            scatter!([newX[1]],[newY[1]],[newZ[1]],label="initial point")
            scatter!([newX[end]],[newY[end]],[newZ[end]],label="end point")
            savefig(string("control/Transfert_",startPoint,"_",endPoint,".png"));

            # Create a 3-row, 1-column layout
            plot_layout = @layout [a; b; c]

            # Create the plots
            p1 = plot(Tvec, newX, title="tangent speed", xlabel="T", ylabel="X", legend=false)
            p2 = plot(Tvec, newY, title="orthogonal speed on ecliptic plane of device", xlabel="T", ylabel="Y", legend=false)
            p3 = plot(Tvec, newZ, title="orthogonal speed", xlabel="T", ylabel="Z", legend=false)

            # Combine the plots into a single window
            plot(p1, p2, p3, layout=plot_layout)
            savefig(string("control/Transfert_2D_",startPoint,"_",endPoint,".png"));

            
            controlVal=[sqrt(sum(sol.control(t).^2))*sol.objective for t in sol.times];
            plot(controlVal,ylim=(0,11));
            savefig(string("control/NORM_Transfert__",startPoint,"_",endPoint,".png"));

            
        end
    end
end
println("finished !!");
CSV.write("result.csv",  Tables.table(objective_output), writeheader=false)
CSV.write("delta_v.csv",  Tables.table(deltas_V), writeheader=false)
CSV.write("mass.csv",  Tables.table(mass), writeheader=false)

mini=0
maxi=0
scatter([0],[0],[0],label="Earth",size=(1920, 1080))
for j in range(1,N-2,step=2)
    states=sols[j,j+2]
    
    X=[x[1] for x in states]
    Y=[x[2] for x in states]
    Z=[x[3] for x in states]
    

    global maxi,mini;

    mini=min(mini,min(minimum(X),minimum(Y),minimum(Z))-0.1)
    maxi=max(maxi,max(maximum(X),maximum(Y),maximum(Z))+0.1)
    plot!(X,Y,Z,label=string(j, " to ",(j+2)),title=string("transfert total ",endPoint),xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi));
end
savefig(string("All_transfert_odd_",startPoint,"_",endPoint,".png"));

mini=0
maxi=0
scatter([0],[0],[0],label="Earth",size=(1920, 1080))
for j in range(2,N-2,step=2)
    states=sols[j,j+2]
    
    X=[x[1] for x in states]
    Y=[x[2] for x in states]
    Z=[x[3] for x in states]
    

    global maxi,mini;

    mini=min(mini,min(minimum(X),minimum(Y),minimum(Z))-0.1)
    maxi=max(maxi,max(maximum(X),maximum(Y),maximum(Z))+0.1)
    plot!(X,Y,Z,label=string(j, " to ",(j+2)),title=string("transfert total ",endPoint),xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi));
end
savefig(string("All_transfert_even_",startPoint,"_",endPoint,".png"));

mini=0
maxi=0
scatter([0],[0],[0],label="Earth",size=(1920, 1080))

gui()
while true
    sleep(0.1)
end
