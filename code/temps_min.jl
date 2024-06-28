println("loading lib");

using OptimalControl
using Plots

println("define pb");

redo=false;

ENDOBS=[-0.167269293387500,
        1.09832700578628,
        0.0450395299341377,
        -3.13998323037230,
        -4.62822855823360,
        -0.450462273954064];

STARTOBS=[0.231357458212836,
          0.922293662828187,
          0.0670351958926985,
          -5.50343680893103,
          -3.61488252858483,
          -0.343544167374632]

Amax=1
muTer=1.0*39.478417604357425
tstart=0.353883781037970
tend=1.409239769931626
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
    ẋ(t) == vcat([vx(t),vy(t),vz(t)],u(t)+gravity([qx(t),qy(t),qz(t)],[0,0,0],muTer))
    #tf → min
    ∫( 0.5sum(u(t).^2) ) → min
end

function F0(X)
    return vcat(X(4:6),gravity(X(1:3),[0,0,0],muTer))
end

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

println("solve pb");
if(redo)
    sol = solve(ocp)
end
println("plot sol");

#plot(sol, size=(1200, 1800))

#gui()

states=[sol.state(x) for x in sol.times]

X=[x[1] for x in states]
Y=[x[2] for x in states]
Z=[x[3] for x in states]

plot(X,Y,Z)
scatter!([STARTOBS[1]],[STARTOBS[2]],[STARTOBS[3]])
scatter!([ENDOBS[1]],[ENDOBS[2]],[ENDOBS[3]])

gui()
while true
    sleep(0.1)
end
