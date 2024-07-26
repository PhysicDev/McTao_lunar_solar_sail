
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

#device related variable :
Isp=3000
Mass=24
FuelMass=33
Thrust=0.02

#ESA smart-1, données à verifier obtenu avec gemini
#Amax ~= 222 µm/s²    6.5 km/s DV

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

fun(t)=t

if(resolve)
    #objective value of the problem (0 mean not computed)
    objective_output=zeros(24,24)
    deltas_V=zeros(24,24)
    
    sols=Array{Any}(undef, 24,24)
    for j in range(1,23,step=1)
        for i in range(j+2,min(j+4,24),step=1)

            global ENDOBS,STARTOBS,startPoint,endPoint,tstart,fun,tend,Tvec,mini,maxi,tspan

            #start and end point of observation
            startPoint=j;
            endPoint=i;
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
                maxia ∈ R, variable
                tf=timef
                ts=time0
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
                ẋ(t) == vcat([vx(t),vy(t),vz(t)],u(t)*maxia+gravity([qx(t),qy(t),qz(t)],[0,0,0],muTer)+ gravity(x(t)[1:3],odeMoon3(t+tstart-time0)[1:3],muLun))
                #tf → min
                maxia → min
            end


            
            #sampling range
            sampling=[50,100,200,400,600];

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

            trueX0(t)=(sat1(t+tstart-time0)*(tend-t+tstart-time0)+sat2(t+tstart-time0)*(t-time0))/(tend-tstart)
            initg = (state=trueX0, control=u0, variable=0)
            #initg=nothing
            println("    ");
            println("solve time pb");#solving problem
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

            println(string("founded solution max acceleration ",sol.objective));

            
            fun(t)=sqrt(sum((sol.control(t)).^2));
            tempDV=integral(fun,time0,timef,10000)*sol.objective;
            println(string("temporary delta V :",tempDV))
            println(string("alternative method : ",sol.objective*(tend-tstart)))




            println("    ")

            println("defining cinetic Pb")
            #defining the problem
            @def ocp begin
                #tf ∈ R, variable
                #maxia = 10
                #useless= 10 , variable

                tf=timef
                ts=time0
                #tend=0.728470473745049
                t ∈ [ ts, tf ], time
                x ∈ R^6, state
                u ∈ R^3, control
                #tf ≥ 0
                sum(u(t).^2)≤10
                #10 ≥ maxia ≥ 0

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
                ẋ(t) == vcat([vx(t),vy(t),vz(t)],u(t)+gravity([qx(t),qy(t),qz(t)],[0,0,0],muTer)+ gravity(x(t)[1:3],odeMoon3(t+tstart-time0)[1:3],muLun))
                #tf → min
                #maxia → min
                ∫( 0.5*(sum(u(t).^2)) ) → min
            end
            cont(t)=sol.control(t)*sol.variable
            initg=(state=sol.state, control=cont)
            println("solving cinetic pb : ")
            sol2=nothing
            #initg=nothing
            
            #sampling range
            sampling=[200,400,600];
            for samp in sampling 
                println(string("grid size : ",samp));
                if(resolve)
                    if(isnothing(initg))#if there are no initial condition do nothing
                        sol2 = solve(ocp,display=false,grid_size=samp)
                    else
                        sol2 = solve(ocp,display=false,grid_size=samp,init=initg)
                    end
                end
                if(0==cmp(sol2.message,"Solve_Succeeded"))#if we solve we set the new intitial solution then we continue
                    println(string("solved  in ",sol2.iterations," iterations || ",sol2.objective))
                    initg=(state=sol.state, control=sol2.control, variable=sol2.variable)
                else#if we don't solve we do nothing
                    println("not solved :(")
                end
            end
            #sol2 = solve(ocp,display=false,grid_size=400,init=initg) #storing solution in another variable in case of convergence failled

            ltpSol=false;            

            fun(t)=sqrt(sum(sol.control(t).^2))
            deltas_V[i,j]=integral(fun,time0,timef,10000)

            if(0==cmp(sol2.message,"Solve_Succeeded"))#if we solve we set use the low fuel solution
                if(deltas_V[i,j]>tempDV)
                    println(string( "DV of low fuel : ",deltas_V[i,j], " greater than DV of low thrust : ",tempDV," ... keeping solution low thrust."))
                    ltpSol=true
                else
                    println(string( "DV of low fuel : ",deltas_V[i,j], " smaller than DV of low thrust : ",tempDV," ... taking solution low fuel."))
                    sol=sol2
                end
                println(string(""))
            else#if we don't solve we keep the low thrust solution
                println("low fuel method failled ... using solution of low thrust problem")
                println(sol2.variable)
                ltpSol=true
            end
            println("   ")

            println(string("plot sol ... optimum : ",sol2.objective));

            #plot(sol, size=(1200, 1800))
            objective_output[i,j]=sol2.objective;#setting objective output variable

            #compute deltaV : 
            if(ltpSol)
                deltas_V[i,j]=deltas_V[i,j]*sol.objective;
            end
            println(string("deltas_V : ",deltas_V[i,j]))

            #now assigning sol to sol2 for plotting and debuging : 
            #sol=sol2

            #deltas_V[i,j]=sol.objective*tend-tstart;#we suppose that we have maximum thrust during the whole transfert
            #gui()

            sols[j,i]=sol;

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

            mini=min(minimum(X),minimum(Y),minimum(Z))-0.1
            maxi=max(maximum(X),maximum(Y),maximum(Z))+0.1


            plot(X,Y,Z,label="trajectory",title=string("transfert from observation ",startPoint," to observation ",endPoint),xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi),size=(1920/2,1080/2));
            
            for t in LowTvec
                St=sol.state(t);
                Ct=sol.control(t)
                if(ltpSol)
                    Ct*=sol.objective;
                end
                plot!([St[1],St[1]+Ct[1]*0.1],[St[2],St[2]+Ct[2]*0.1],[St[3],St[3]+Ct[3]*0.1],label=false,color=:red)
            end
            scatter!([STARTOBS[1]],[STARTOBS[2]],[STARTOBS[3]],label="starting point")
            scatter!([ENDOBS[1]],[ENDOBS[2]],[ENDOBS[3]],label="final point")
            scatter!([0],[0],[0],label="Earth")


            #ploting start and end orbit:
            plot!(sat1,vars=(1,2,3),label="starting obrit",color=:green)
            plot!(sat2,vars=(1,2,3),label="ending obrit",color=:orange)
            



            #title("transfert from observation "+startPoint+" to observation "+endPoint)
            savefig(string("images/Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));

            plot(sol, size=(1200, 1800))
            savefig(string("solGraph/Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));

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
            savefig(string("control/Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));

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
            savefig(string("control/NORM_Transfert_Amax_",Amax,"_",startPoint,"_",endPoint,".png"));

            
        end
    end  
end
println("finished !!");
CSV.write("result.csv",  Tables.table(objective_output), writeheader=false)
CSV.write("delta_v.csv",  Tables.table(deltas_V), writeheader=false)

mini=0
maxi=0
scatter([0],[0],[0],label="Earth",size=(1920, 1080))
for j in range(1,21,step=2)
    states=sols[j,j+2]
    
    X=[x[1] for x in states]
    Y=[x[2] for x in states]
    Z=[x[3] for x in states]
    

    global maxi,mini;

    mini=min(mini,min(minimum(X),minimum(Y),minimum(Z))-0.1)
    maxi=max(maxi,max(maximum(X),maximum(Y),maximum(Z))+0.1)
    plot!(X,Y,Z,label=string(j, " to ",(j+2)),title=string("transfert total ",endPoint),xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi));
end
savefig(string("All_transfert_odd_",Amax,"_",startPoint,"_",endPoint,".png"));

mini=0
maxi=0
scatter([0],[0],[0],label="Earth",size=(1920, 1080))
for j in range(2,22,step=2)
    states=sols[j,j+2]
    
    X=[x[1] for x in states]
    Y=[x[2] for x in states]
    Z=[x[3] for x in states]
    

    global maxi,mini;

    mini=min(mini,min(minimum(X),minimum(Y),minimum(Z))-0.1)
    maxi=max(maxi,max(maximum(X),maximum(Y),maximum(Z))+0.1)
    plot!(X,Y,Z,label=string(j, " to ",(j+2)),title=string("transfert total ",endPoint),xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi));
end
savefig(string("All_transfert_even_",Amax,"_",startPoint,"_",endPoint,".png"));

mini=0
maxi=0
scatter([0],[0],[0],label="Earth",size=(1920, 1080))
for j in range(1,19,step=3)
    states=sols[j,j+3]
    
    X=[x[1] for x in states]
    Y=[x[2] for x in states]
    Z=[x[3] for x in states]
    

    global maxi,mini;

    mini=min(mini,min(minimum(X),minimum(Y),minimum(Z))-0.1)
    maxi=max(maxi,max(maximum(X),maximum(Y),maximum(Z))+0.1)
    plot!(X,Y,Z,label=string(j, " to ",(j+2)),title=string("transfert total ",endPoint),xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi));
end
savefig(string("All_transfert_1_on_3_",Amax,"_",startPoint,"_",endPoint,".png"));

mini=0
maxi=0
scatter([0],[0],[0],label="Earth",size=(1920, 1080))
for j in range(1,17,step=4)
    states=sols[j,j+4]
    
    X=[x[1] for x in states]
    Y=[x[2] for x in states]
    Z=[x[3] for x in states]
    

    global maxi,mini;

    mini=min(mini,min(minimum(X),minimum(Y),minimum(Z))-0.1)
    maxi=max(maxi,max(maximum(X),maximum(Y),maximum(Z))+0.1)
    plot!(X,Y,Z,label=string(j, " to ",(j+2)),title=string("transfert total ",endPoint),xlim=(mini, maxi),ylim=(mini, maxi),zlim=(mini, maxi));
end
savefig(string("All_transfert_1_on_4.png"));

gui()
while true
    sleep(0.1)
end
