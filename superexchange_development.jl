
include("constants.jl")

#driving_force = parse(Float64,ARGS[1])
the_parameter = driving_force

#is12 = parse(Float64,ARGS[1])
#the_parameter = is12
#is23=is12

#tunneling_energy = parse(Float64,ARGS[1])
#the_parameter = tunneling_energy

function esys(x::Float64,hsys::Array{Float64,2},eval::Array{Float64,1},evec::Array{Float64,2})
    hsys[1,1]=(0.5*mass*omega^(2.0))*x^(2.0)+subotnik*x
    hsys[2,2]=(0.5*mass*omega^(2.0))*(x+driving_force/(2.0*subotnik))^(2.0)+tunneling_energy
    hsys[3,3]=(0.5*mass*omega^(2.0))*x^(2.0)-subotnik*x-driving_force
    hsys[1,2]=hsys[2,1]=is12
    hsys[1,3]=hsys[3,1]=0.00
    hsys[2,3]=hsys[3,2]=is23
    eval[:],evec[:,:]=eig(hsys)
    index=sortperm(eval)
    eval[:]=eval[index]
    evec[:,:]=evec[:,index]
    nothing    
end

function force_matrix(x::Float64,fsys::Array{Float64,2})
    fsys[1,1]=(mass*omega^(2.0))*x+subotnik
    fsys[2,2]=(mass*omega^(2.0))*(x+driving_force/(2.0*subotnik))
    fsys[3,3]=(mass*omega^(2.0))*x-subotnik
    fsys[1,2]=fsys[2,1]=0.0
    fsys[1,3]=fsys[3,1]=0.0
    fsys[2,3]=fsys[3,2]=0.0
    nothing
end 

function elemental_force(f::Array{Float64,1},bra::Array{Float64,1},ket::Array{Float64,1},fsys::Array{Float64,2})
    f[1]=-dot(bra,fsys*ket)
    nothing
end

function ehrenfest_decohere_force(v::Float64,f::Array{Float64,1},eval::Array{Float64,1},nf::Array{Float64,2},rho::Array{Complex{Float64},2},dm::Array{Complex{Float64},2})
    f[1]=real(sum(rho.*nf)-(mass/v)*sum(diag(dm).*eval))
    nothing
end

function force(ftype::Int64,cstate::Int64,v::Float64,f::Array{Float64,1},eval::Array{Float64,1},evec::Array{Float64,2},fsys::Array{Float64,2},nf::Array{Float64,2},rho::Array{Complex{Float64},2},dm::Array{Complex{Float64},2})
    if ftype==1
        elemental_force(f,evec[:,cstate],evec[:,cstate],fsys)
    end
    if ftype==2
        ehrenfest_decohere_force(v,f,eval,nf,rho,dm)
    end
    nothing
end

function classical(ftype::Int64,cstate::Int64,xv::Array{Float64,1},xvc::Array{Float64,1},dx::Array{Float64,1},dv::Array{Float64,1},f::Array{Float64,1},eval::Array{Float64,1},hsys::Array{Float64,2},evec::Array{Float64,2},fsys::Array{Float64,2},nf::Array{Float64,2},rho::Array{Complex{Float64},2},dm::Array{Complex{Float64},2})
    xvc[1]=xv[1]
    xvc[2]=xv[2]
    ## rk4:1 ##
    esys(xv[1],hsys,eval,evec)
    force_matrix(xv[1],fsys)
    force(ftype,cstate,xv[2],f,eval,evec,fsys,nf,rho,dm)
    dv[1]=(dt/mass)*(f[1]-mass*friction*xv[2])
    dx[1]=dt*xv[2]
    xv[1]=xvc[1]+dx[1]/2.0
    xv[2]=xvc[2]+dv[1]/2.0
    ## rk4:2 ##
    esys(xv[1],hsys,eval,evec)
    force_matrix(xv[1],fsys)
    force(ftype,cstate,xv[2],f,eval,evec,fsys,nf,rho,dm)
    dv[2]=(dt/mass)*(f[1]-mass*friction*xv[2])
    dx[2]=dt*xv[2]
    xv[1]=xvc[1]+dx[2]/2.0
    xv[2]=xvc[2]+dv[2]/2.0
    ## rk4:3 ##
    esys(xv[1],hsys,eval,evec)
    force_matrix(xv[1],fsys)
    force(ftype,cstate,xv[2],f,eval,evec,fsys,nf,rho,dm)
    dv[3]=(dt/mass)*(f[1]-mass*friction*xv[2])
    dx[3]=dt*xv[2]
    xv[1]=xvc[1]+dx[3]
    xv[2]=xvc[2]+dv[3]
    ## rk4:4 ##
    esys(xv[1],hsys,eval,evec)
    force_matrix(xv[1],fsys)
    force(ftype,cstate,xv[2],f,eval,evec,fsys,nf,rho,dm)
    dv[4]=(dt/mass)*(f[1]-mass*friction*xv[2])
    dx[4]=dt*xv[2]
    ### end integration ###
    xv[1]=xvc[1]+(dx[1]+2.0*(dx[2]+dx[3])+dx[4])/6.0
    xv[2]=xvc[2]+(dv[1]+2.0*(dv[2]+dv[3])+dv[4])/6.0
    nothing
end

function thermostat(ftype::Int64,cstate::Int64,xv::Array{Float64,1},r::Array{Float64,1},z::Array{Float64,1},t::Array{Float64,1},f::Array{Float64,1},eval::Array{Float64,1},hsys::Array{Float64,2},evec::Array{Float64,2},fsys::Array{Float64,2},nf::Array{Float64,2},rho::Array{Complex{Float64},2},dm::Array{Complex{Float64},2})
    r[:]=randn(4)*langevin_sigma
    #r[:]=[0.001 -0.001 0.001 -0.001]*langevin_sigma
    z[1]=r[1]
    z[2]=dt*(r[1]/2.0+r[2]/(2.0*sqrt(3.0)))
    z[3]=dt^(2.0)*(r[1]/6.0+(sqrt(3.0)/12.0)*r[2]+r[3]/(12.0*sqrt(5.0)))
    z[4]=dt^(3.0)*(r[1]/24.0+(sqrt(3.0)/40.0)*r[2]+r[3]/(24.0*sqrt(5.0))+r[4]/(120.0*sqrt(7.0)))
    esys(xv[1]+2.0*dx_num_deriv,hsys,eval,evec)
    force_matrix(xv[1]+2.0*dx_num_deriv,fsys)
    force(ftype,cstate,xv[2],f,eval,evec,fsys,nf,rho,dm)
    t[1]=-f[1]
    esys(xv[1]+dx_num_deriv,hsys,eval,evec)
    force_matrix(xv[1]+dx_num_deriv,fsys)
    force(ftype,cstate,xv[2],f,eval,evec,fsys,nf,rho,dm)
    t[2]=8.0*f[1]
    esys(xv[1]-dx_num_deriv,hsys,eval,evec)
    force_matrix(xv[1]-dx_num_deriv,fsys)
    force(ftype,cstate,xv[2],f,eval,evec,fsys,nf,rho,dm)
    t[3]=-8.0*f[1]
    esys(xv[1]-2.0*dx_num_deriv,hsys,eval,evec)
    force_matrix(xv[1]-2.0*dx_num_deriv,fsys)
    force(ftype,cstate,xv[2],f,eval,evec,fsys,nf,rho,dm)
    t[4]=f[1]
    hessian=-(t[1]+t[2]+t[3]+t[4])/(12.0*dx_num_deriv)
    ### position and velocity ###
    xv[1]=xv[1]+z[2]-friction*z[3]+(-hessian/mass+friction^(2.0))*z[4]
    xv[2]=xv[2]+z[1]-friction*z[2]+(-hessian/mass+friction^(2.0))*z[3]+(2.0*friction*hessian/mass-friction^(3.0))*z[4]
    nothing
end

function norm(evec::Array{Float64,2},evecc::Array{Float64,2},evecn::Array{Float64,2})
    evecn[:,:]=evec
    for i=1:nstates
        if sign(dot(evecc[:,i],evecn[:,i]))==-1.0
            evecn[:,i]=-evecn[:,i]
        end
    end
    evec[:,:]=evecn
    evecc[:,:]=evecn
    nothing
end

function nac(f::Array{Float64,1},d::Array{Float64,2},nf::Array{Float64,2},df::Array{Float64,2},cstate::Int64,evec::Array{Float64,2},eval::Array{Float64,1},fsys::Array{Float64,2})
    for i=1:nstates
        for j=1:i
            elemental_force(f,evec[:,i],evec[:,j],fsys)
            nf[i,j]=f[1]
            nf[j,i]=nf[i,j]
            if i==j
                d[i,i]=0.0
            else
                d[i,j]=-nf[i,j]/(eval[j]-eval[i])
                d[j,i]=-d[i,j]
            end
        end
    end
    df[:,:]=nf-nf[cstate,cstate]*eye(nstates)
    nothing
end

function moments(xm::Array{Complex{Float64},2},pm::Array{Complex{Float64},2},v::Float64,cstate::Int64,d::Array{Float64,2},hd::Array{Float64,2},df::Array{Float64,2},dmat::Array{Complex{Float64},2})
    xm[:,:]=-1.0im*(hd*xm-xm*hd)+pm/mass-v*(d*xm-xm*d)
    pm[:,:]=-1.0im*(hd*pm-pm*hd)+(df*dmat+dmat*df)/2.0-v*(d*pm-pm*d)
    xm[:,:]=dt*(xm-xm[cstate,cstate]*eye(nstates))
    pm[:,:]=dt*(pm-pm[cstate,cstate]*eye(nstates))
    nothing
end

function integrating_moments(xm::Array{Complex{Float64},2},pm::Array{Complex{Float64},2},xmc::Array{Complex{Float64},2},pmc::Array{Complex{Float64},2},v::Float64,cstate::Int64,d::Array{Float64,2},hd::Array{Float64,2},df::Array{Float64,2},dmat::Array{Complex{Float64},2},ddp1::Array{Complex{Float64},2},ddp2::Array{Complex{Float64},2},ddp3::Array{Complex{Float64},2},ddp4::Array{Complex{Float64},2},ddx1::Array{Complex{Float64},2},ddx2::Array{Complex{Float64},2},ddx3::Array{Complex{Float64},2},ddx4::Array{Complex{Float64},2})
    xmc[:,:]=xm
    pmc[:,:]=pm
    ## rk4:1 ##
    moments(xm,pm,v,cstate,d,hd,df,dmat)
    ddp1[:,:]=pm
    ddx1[:,:]=xm
    xm[:,:]=xmc+ddx1/2.0
    pm[:,:]=pmc+ddp1/2.0
    ## rk4:2 ##
    moments(xm,pm,v,cstate,d,hd,df,dmat)
    ddp2[:,:]=pm
    ddx2[:,:]=xm
    xm[:,:]=xmc+ddx2/2.0
    pm[:,:]=pmc+ddp2/2.0
    ## rk4:3 ##
    moments(xm,pm,v,cstate,d,hd,df,dmat)
    ddp3[:,:]=pm
    ddx3[:,:]=xm
    xm[:,:]=xmc+ddx3
    pm[:,:]=pmc+ddp3
    ## rk4:4 ##
    moments(xm,pm,v,cstate,d,hd,df,dmat)
    ddp4[:,:]=pm
    ddx4[:,:]=xm
    ### end integration ###
    xm[:,:]=xmc+(ddx1+2.0*(ddx2+ddx3)+ddx4)/6.0
    pm[:,:]=pmc+(ddp1+2.0*(ddp2+ddp3)+ddp4)/6.0
    nothing
end

function quantum(c::Array{Complex{Float64},1},h::Array{Complex{Float64},2},dmat::Array{Complex{Float64},2},dc1::Array{Complex{Float64},1},dc2::Array{Complex{Float64},1},dc3::Array{Complex{Float64},1},dc4::Array{Complex{Float64},1})
    dc1[:]=-1.0im*dt*h*c
    dc2[:]=-1.0im*dt*h*(c+dc1/2.0)
    dc3[:]=-1.0im*dt*h*(c+dc2/2.0)
    dc4[:]=-1.0im*dt*h*(c+dc3)
    c[:]=c+(dc1+2.0*(dc2+dc3)+dc4)/6.0
    ### end integration ###
    dmat[:,:]=c*ctranspose(c)
    nothing
end

function quantum_dm(dmat::Array{Complex{Float64},2},h::Array{Complex{Float64},2},ddmat1::Array{Complex{Float64},2},ddmat2::Array{Complex{Float64},2},ddmat3::Array{Complex{Float64},2},ddmat4::Array{Complex{Float64},2})
    ddmat1[:,:]=-1.0im*dt*(h*dmat-dmat*h)
    ddmat2[:,:]=-1.0im*dt*(h*(dmat+ddmat1/2.0)-(dmat+ddmat1/2.0)*h)
    ddmat3[:,:]=-1.0im*dt*(h*(dmat+ddmat2/2.0)-(dmat+ddmat2/2.0)*h)
    ddmat4[:,:]=-1.0im*dt*(h*(dmat+ddmat3)-(dmat+ddmat3)*h)
    ### end integration ###
    dmat[:,:]=dmat+(ddmat1+2.0*(ddmat2+ddmat3)+ddmat4)/6.0
    nothing
end

function gfsh(cnastate::Array{Int64,1},xv::Array{Float64,1},deltap::Array{Float64,1},pbefore::Array{Float64,1},pafter::Array{Float64,1},positives::Array{Float64,1},hopping::Array{Float64,1},eval::Array{Float64,1},xm::Array{Complex{Float64},2},pm::Array{Complex{Float64},2},dmat::Array{Complex{Float64},2})
    pafter[:]=real(diag(dmat))
    deltap[:]=pafter-pbefore
    deltap[isnan.(deltap)]=0.0
    if deltap[cnastate[1]]>0.0
        cnastate[2]=cnastate[1]
        cnastate[3]=cnastate[1]
    else
        positives[:]=deltap
        positives[positives.<0.0]=0.0
        hopping[:]=cumsum(-deltap[cnastate[1]]*positives/(pafter[cnastate[1]]*sum(positives)))
        hopping[isnan.(hopping)]=0.0
        index=searchsortedfirst(hopping,rand())
        #index=searchsortedfirst(hopping,0.2)
        if index==nstates+1
            cnastate[2]=cnastate[1]
            cnastate[3]=cnastate[1]
        else
            cnastate[2]=states[index]
            cnastate[3]=states[index]
            if mass*abs2(xv[2])/2.0+eval[cnastate[1]]<eval[cnastate[2]]
                cnastate[2]=cnastate[1]
                xv[2]=-xv[2]
            else
                xv[2]=sign(xv[2])*sqrt(abs2(xv[2])-(2.0/mass)*(eval[cnastate[2]]-eval[cnastate[1]]))
                xm[:,:]=0.0
                pm[:,:]=0.0
            end
        end
    end
    pbefore[:]=pafter
    nothing
end

function decohere_time(inv_tau::Array{Float64,1},i::Int64,j::Int64,v::Float64,eval::Array{Float64,1})
    inv_tau[1]=abs(eval[i]-eval[j])*(c_truhlar+(mass*abs2(v))/(2.0*e0_truhlar))^(-1.0)
    nothing
end

function demixing_matrix(cstate::Int64,v::Float64,inv_tau::Array{Float64,1},eval::Array{Float64,1},rho::Array{Complex{Float64},2},h::Array{Complex{Float64},2},cm::Array{Complex{Float64},2},dm::Array{Complex{Float64},2})
    cm[:,:]=-1.0im*(h*rho-rho*h)
    dm[:,:]=zeros(Complex{Float64},3,3)
    st=0.0
    for k in states[states.!=cstate]
        decohere_time(inv_tau,cstate,k,v,eval)
        st+=real(rho[k,k]*inv_tau[1])
    end
    for i=1:nstates
        for j=1:i
            if i==j
                if i==cstate
                    dm[i,i]=0.5*st
                else
                    decohere_time(inv_tau,i,cstate,v,eval)
                    dm[i,i]=-0.5*rho[i,i]*inv_tau[1]
                end
            else
                if i!=cstate && j!=cstate
                    decohere_time(inv_tau,i,cstate,v,eval)
                    dm[i,j]=inv_tau[1]
                    decohere_time(inv_tau,j,cstate,v,eval)
                    dm[i,j]=-0.5*(dm[i,j]+inv_tau[1])*rho[i,j]
                end
                if i==cstate && j!=cstate
                    decohere_time(inv_tau,j,cstate,v,eval)
                    dm[i,j]=0.5*((1.0/rho[cstate,cstate])*st-inv_tau[1])*rho[i,j]
                end
                if i!=cstate && j==cstate
                    decohere_time(inv_tau,i,cstate,v,eval)
                    dm[i,j]=0.5*((1.0/rho[cstate,cstate])*st-inv_tau[1])*rho[i,j]
                end
            end
        end
    end
    dm[:,:]=dm+ctranspose(dm)
    rho[:,:]=dt*(cm+dm)
    nothing
end

function integrating_demixing_matrix(cstate::Int64,v::Float64,inv_tau::Array{Float64,1},eval::Array{Float64,1},rho::Array{Complex{Float64},2},rhoc::Array{Complex{Float64},2},h::Array{Complex{Float64},2},cm::Array{Complex{Float64},2},dm::Array{Complex{Float64},2},drho1::Array{Complex{Float64},2},drho2::Array{Complex{Float64},2},drho3::Array{Complex{Float64},2},drho4::Array{Complex{Float64},2})
    rhoc[:,:]=rho
    ## rk4:1 ##
    demixing_matrix(cstate,v,inv_tau,eval,rho,h,cm,dm)
    drho1[:,:]=rho
    rho[:,:]=rhoc+drho1/2.0
    ## rk4:2 ##
    demixing_matrix(cstate,v,inv_tau,eval,rho,h,cm,dm)
    drho2[:,:]=rho
    rho[:,:]=rhoc+drho2/2.0
    ## rk4:3 ##
    demixing_matrix(cstate,v,inv_tau,eval,rho,h,cm,dm)
    drho3[:,:]=rho
    rho[:,:]=rhoc+drho3
    ## rk4:4 ##
    demixing_matrix(cstate,v,inv_tau,eval,rho,h,cm,dm)
    drho4[:,:]=rho
    ### end integration ###
    rho[:,:]=rhoc+(drho1+2.0*(drho2+drho3)+drho4)/6.0
    nothing
end

function decohere_augmented_algorithm(cstate::Int64,rates::Array{Float64,1},nf::Array{Float64,2},df::Array{Float64,2},c::Array{Complex{Float64},1},xm::Array{Complex{Float64},2},pm::Array{Complex{Float64},2})
    for pstate in states[states.!=cstate]
        rates[1]=real(-(dt/2.0)*df[pstate,pstate]*xm[pstate,pstate])
        rates[2]=real(-rates[1]-(2.0*dt)*abs(nf[cstate,pstate]*xm[pstate,pstate]))
        rd=rand()
        #rd=0.2
        if rd<rates[2]
            c[:]=c/sqrt(1.0-abs2(c[pstate]))
            c[pstate]=0.0
            xm[pstate,:]=xm[:,pstate]=0.0
            pm[pstate,:]=pm[:,pstate]=0.0
        else
            if rd<rates[1]
                xm[pstate,:]=xm[:,pstate]=0.0
                pm[pstate,:]=pm[:,pstate]=0.0
            end
        end
    end
    nothing
end

function decohere_simple_collapse_on_successful_hop(cnastate::Array{Int64,1},c::Array{Complex{Float64},1})
    if cnastate[2]!=cnastate[1]
        c[:]=zeros(nstates)
        c[cnastate[2]]=1.0
    end
    nothing
end
    
function decohere_simple_collapse_on_attempted_hop(cnastate::Array{Int64,1},c::Array{Complex{Float64},1})
    if cnastate[3]!=cnastate[1]
        c[:]=zeros(nstates)
        c[cnastate[3]]=1.0
    end
    nothing
end

function main()
    
    ### allocate memory ###
    xv=Array{Float64}(2)
    xvc=Array{Float64}(2)
    hsys=Array{Float64}((nstates,nstates))
    eval=Array{Float64}(nstates)
    evec=Array{Float64}((nstates,nstates))
    fsys=Array{Float64}((nstates,nstates))
    f=Array{Float64}(1)
    dx=Array{Float64}(4)
    dv=Array{Float64}(4)
    r=Array{Float64}(4)
    z=Array{Float64}(4)
    t=Array{Float64}(4)
    evecc=Array{Float64}((nstates,nstates))
    evecn=Array{Float64}((nstates,nstates))
    d=Array{Float64}((nstates,nstates))
    nf=Array{Float64}((nstates,nstates))
    h=Array{Complex{Float64}}((nstates,nstates))
    c=Array{Complex{Float64}}(nstates)
    dc1=Array{Complex{Float64}}(nstates)
    dc2=Array{Complex{Float64}}(nstates)
    dc3=Array{Complex{Float64}}(nstates)
    dc4=Array{Complex{Float64}}(nstates)
    dmat=Array{Complex{Float64}}((nstates,nstates))
    cnastate=Array{Int64}(3)
    pbefore=Array{Float64}(nstates)
    pafter=Array{Float64}(nstates)
    deltap=Array{Float64}(nstates)
    positives=Array{Float64}(nstates)
    hopping=Array{Float64}(nstates)
    rates=Array{Float64}(2)
    diabatp=Array{Float64}(nstates)
    timestep::Int64 = 0
    tsindex::Int64 = 1
    answer=Array{Float64}((nstates+1,anas))
    
    ### truhlar decay-of-mixing density of matrix ###
    rho=Array{Complex{Float64}}((nstates,nstates))
    rhoc=Array{Complex{Float64}}((nstates,nstates))
    drho1=Array{Complex{Float64}}((nstates,nstates))
    drho2=Array{Complex{Float64}}((nstates,nstates))
    drho3=Array{Complex{Float64}}((nstates,nstates))
    drho4=Array{Complex{Float64}}((nstates,nstates))
    dm=Array{Complex{Float64}}((nstates,nstates))
    cm=Array{Complex{Float64}}((nstates,nstates))
    inv_tau=Array{Float64}(1)
    
    ### liouville equation to integrate density of matrix ###
    ddmat1=Array{Complex{Float64}}((nstates,nstates))
    ddmat2=Array{Complex{Float64}}((nstates,nstates))
    ddmat3=Array{Complex{Float64}}((nstates,nstates))
    ddmat4=Array{Complex{Float64}}((nstates,nstates)) 
    
    ### subotnik decoherence ###
    df=Array{Float64}((nstates,nstates))
    xm=Array{Complex{Float64}}((nstates,nstates))
    pm=Array{Complex{Float64}}((nstates,nstates))
    xmc=Array{Complex{Float64}}((nstates,nstates))
    pmc=Array{Complex{Float64}}((nstates,nstates))
    hd=Array{Float64}((nstates,nstates))
    ddx1=Array{Complex{Float64}}((nstates,nstates))
    ddx2=Array{Complex{Float64}}((nstates,nstates))
    ddx3=Array{Complex{Float64}}((nstates,nstates))
    ddx4=Array{Complex{Float64}}((nstates,nstates))
    ddp1=Array{Complex{Float64}}((nstates,nstates))
    ddp2=Array{Complex{Float64}}((nstates,nstates))
    ddp3=Array{Complex{Float64}}((nstates,nstates))
    ddp4=Array{Complex{Float64}}((nstates,nstates))
    
    ### initialize trajectory ###
    xv[1]=mean_x+randn()*sigma_x
    xv[2]=(mean_p+randn()*sigma_p)/mass
    #xv[1]=mean_x+0.1*sigma_x
    #xv[2]=(mean_p+0.1*sigma_p)/mass
    esys(xv[1],hsys,eval,evec)
    evecc[:,:]=evec
    cnastate[1]=searchsortedfirst(cumsum(abs2.(evecc[1,:])),rand())
    #cnastate[1]=searchsortedfirst(cumsum(abs2.(evecc[1,:])),0.1)
    c[:]=evecc[1,:]
    dmat[:,:]=c*ctranspose(c)
    pbefore[:]=real(diag(dmat))
    rho[:,:]=dmat
    
    ### begin trajectory ###
    for timestep=1:tots
        classical(ftype,cnastate[1],xv,xvc,dx,dv,f,eval,hsys,evec,fsys,nf,rho,dm)
        thermostat(ftype,cnastate[1],xv,r,z,t,f,eval,hsys,evec,fsys,nf,rho,dm)
        esys(xv[1],hsys,eval,evec)
        norm(evec,evecc,evecn)
        force_matrix(xv[1],fsys)
        nac(f,d,nf,df,cnastate[1],evec,eval,fsys)
        hd[:,:]=diagm(eval)
        #integrating_moments(xm,pm,xmc,pmc,xv[2],cnastate[1],d,hd,df,dmat,ddp1,ddp2,ddp3,ddp4,ddx1,ddx2,ddx3,ddx4)
        h[:,:]=hd-1.0*im*xv[2]*d
        #quantum(c,h,dmat,dc1,dc2,dc3,dc4)
        #quantum_dm(dmat,h,ddmat1,ddmat2,ddmat3,ddmat4)
        integrating_demixing_matrix(cnastate[1],xv[2],inv_tau,eval,rho,rhoc,h,cm,dm,drho1,drho2,drho3,drho4)
        gfsh(cnastate,xv,deltap,pbefore,pafter,positives,hopping,eval,xm,pm,rho) # dmat is coherent, rho is decoherent
        #decohere_simple_collapse_on_successful_hop(cnastate,c)
        #decohere_augmented_algorithm(cnastate[1],rates,nf,df,c,xm,pm)
        diabatp[:]=abs2.(evec)[:,cnastate[2]]
        if timestep in print_indices
            @printf "%.02f %.04f %.04f %.04f %.04f %.04f %.04f %.04f %.04f\n" timestep*dt xv[1] eval[cnastate[2]] abs(rho[1,1]) abs(rho[2,2]) abs(rho[3,3]) abs(rho[1,2]) abs(rho[1,3]) abs(rho[2,3])
        end
        #if timestep in print_indices
        #    @printf "%.02f %.04f %.04f %.04f %.04f %.04f %.04f %.04f %.04f\n" timestep*dt xv[1] eval[cnastate[2]] abs(dmat[1,1]) abs(dmat[2,2]) abs(dmat[3,3]) abs(dmat[1,2]) abs(dmat[1,3]) abs(dmat[2,3])
        #end
        #if timestep in tsindices
        #    answer[1,tsindex]=timestep*dt
        #    answer[2:(nstates+1),tsindex]=diabatp
        #    tsindex += 1
        #end
        cnastate[1]=cnastate[2]
        nothing
    end
    #@printf "%.f\n" the_parameter
    #for i=1:anas
    #    @printf "%.02f %.04f %.04f %.04f\n" answer[1,i] answer[2,i] answer[3,i] answer[4,i]
    #end
end

main()
