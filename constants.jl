#!/usr/bin/env julia

### timesteps ###
const dt=1.25
const tots=10^(7) # 7
const anas=10^(2) # 2
const tsindices=range(1,convert(Int,tots/anas),anas)
const print_indices=range(1,convert(Int,tots/10^(4)),10^(4))

### states ###
const nstates=3
const states=range(1,nstates)

### model parameters ###
const ifactor=1.000
const mass=1.000
const omega=4.375*10.0^(-5.0)
const driving_force=2.400*10.0^(-2.0)
const tunneling_energy=10.0*10.0^(-2.0)
const reorganization=2.390*10.0^(-2.0)
const subotnik=sqrt(reorganization*mass*omega^2.0/2.0)
const is12=sqrt(2.500*10.0^(-6.0))
const is23=ifactor*is12

### force type ###
const ftype=1 # 1 = single surface, 2 = Ehrenfest force with decoherent force

### truhlar decoherence parameters ###
const e0_truhlar=0.1
const c_truhlar=1.0 # changing from 10.0

### langevin parameters ###
const temp=9.500*10.0^(-4.0)
const friction=1.500*10.0^(-4.0)
const langevin_sigma=sqrt(2.000*friction*temp*dt/mass)
const dx_num_deriv=1.000*10.0^(-4.0)

### initial classical-coordinate conditions ###
const sigma_x=sqrt(temp/(mass*omega^2.0))
const sigma_p=sqrt(mass*temp)
const mean_x=-subotnik/(mass*omega^(2.0))
const mean_p=0.0
