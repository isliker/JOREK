!> Push relativistic guiding centre(s) in JOREK fields
!>
!> Compile with `make ex7_jorek`
!> Run with `./ex7_jorek < JOREK_namelist`

program ex7_jorek

use particle_tracer
use mod_particle_io
use mod_particle_diagnostics
use mod_fields_linear   
use mod_fields_hermite_birkhoff
use mod_gc_relativistic

use mod_kinetic_relativistic

use hdf5_io_module
                  
implicit none

! Set up the simulation variables
real(kind=8)                      :: timesteps(1) = [1.25d-11] ![3.5723d-13] ! 
real(kind=8)                      :: target_time, t
real*8                            :: energy !< Initial energy in eV (including rest energy)
real*8                            :: ksi    !< Initial cosine of pitch-angle
real*8                            :: chi
integer(kind=4)                   :: n_part, i, j, k, n_steps, ifail, n_lost, i_gc
logical                           :: restart
type(diag_print_kinetic_energy)   :: print_kinetic_energy
type(write_particle_diagnostics)  :: diag
    ! type `write_particle_diagnostics` is defined in `mod_particle_diagnostics`,
    ! extends the type `io_action` and adds array `only`,
    ! `io_action` is defined in `mod_io_actions` and extends type `action`,
    ! `action` is defined in `mod_event`
type(particle_gc_relativistic)    :: particle_in, particle_out
type(particle_kinetic_relativistic)    :: particle_in_L, particle_out_L
type(particle_gc)                 :: particle_out_gc
real(kind=8)                      :: psi, U, gyro_angle
real(kind=8),dimension(3)         :: E, B
! For interp_PRZ
!real*8 :: P(1), P_s(1), P_t(1), P_phi(1), P_time(1)
real*8 :: R, R_s, R_t, Z, Z_s, Z_t

real*8 :: duration, diag_step, anini
integer(hid_t) :: h5_file_id
integer(kind=4) :: ierr, nini, ini
real*8, dimension(:), allocatable :: xini, yini, zini, cosini, Etotalini 
real*8, dimension(:), allocatable :: Rini, phiini

i_gc = 0 ! guiding center (1) or Lorentz (0)

call sim%initialize(num_groups=1)

! Set up the diagnostics output
diag = write_particle_diagnostics(filename='diag.h5') !,only=[1,2,6,12,13,14,15]) ! store total and kinetic energies, p_phi, ielm, phi, R, Z

restart = .false. !

if (.not. restart) then
  call hdf5_open('ini_cond.h5', h5_file_id, ierr)
  call hdf5_real_reading(   h5_file_id, anini,      '/position/nini')
  nini = nint (anini)
  allocate (xini(nini), yini(nini), zini(nini))
  allocate (Rini(nini), phiini(nini))
  allocate (cosini(nini), Etotalini(nini) )
  call hdf5_array1d_reading(h5_file_id, xini,      '/position/xini')
  call hdf5_array1d_reading(h5_file_id, yini,      '/position/yini')
  call hdf5_array1d_reading(h5_file_id, zini,      '/position/zini')
  call hdf5_array1d_reading(h5_file_id, Rini,      '/position/Rini')
  call hdf5_array1d_reading(h5_file_id, phiini,    '/position/phiini_JOREK')
  call hdf5_array1d_reading(h5_file_id, cosini,    '/velocity/cosini')
  call hdf5_array1d_reading(h5_file_id, Etotalini, '/velocity/Etotalini_eV')
  call hdf5_close(h5_file_id)
  write(6,*) 'xini', xini(1), xini(nini)
  write(6,*) 'yini', yini(1), yini(nini)
  write(6,*) 'zini', zini(1), zini(nini)
  write(6,*) 'Rini', Rini(1), Rini(nini)
  write(6,*) 'phiini', phiini(1), phiini(nini)
  write(6,*) 'cosini', cosini(1), cosini(nini)
  write(6,*) 'Etotalini', Etotalini(1), Etotalini(nini)

  sim%time = 0.03251098d0 ! 2.5d-3 ! 1.d-7 !   ! start time in sec

  ! Allocate a group and particle(s) of type particle_gc_relativistic
  n_part = 48
  !allocate(particle_gc_relativistic::sim%groups(1)%particles(n_part))
  sim%groups(1)%mass = 5.4857990907016d-4 !< particle mass in AMU: electron mass
  allocate(particle_kinetic_relativistic::sim%groups(1)%particles(n_part))
  ! Set events to write output data and stop the simulation.
  ! One can use read_jorek_fields_interp_linear or read_jorek_fields_interp_hermite_birkhoff,
  ! and i=-1 (to read jorek_restart.h5 and keep this field at all time) or i=last_file_before_time(sim%time)
  ! (to read a sequel of jorekXXXXX.h5 files and use time-evolving fields)
  duration = 7.0d-4
  !duration = 0.056d-3
  !duration = duration/10.0d0
  diag_step = duration/100.0d0
  !events = [event(read_jorek_fields_interp_linear(i=last_file_before_time(sim%time))), & 
  !          event(diag, start = sim%time, step = diag_step),         &
  !          event(stop_action(), start = sim%time + duration)]
  events = [event(read_jorek_fields_interp_hermite_birkhoff(i=last_file_before_time(sim%time))), & 
            event(diag, start = sim%time, step = diag_step),         &
            event(stop_action(), start = sim%time + duration)]

  ! Run first event to read the JOREK fields
  call with(sim, events, at=0.d0)

  do ini = 1, n_part
     select type (p=>sim%groups(1)%particles(ini))
     type is (particle_kinetic_relativistic)
        particle_in%q = -1
        !p%x = [1.90d0, 0.d0, 0.d0] ! = (R, Z, phi) in [m, m, rad]
        particle_in%x = [Rini(ini), zini(ini), phiini(ini)] ! = (R, Z, phi) in [m, m, rad]
        call find_RZ(sim%fields%node_list, sim%fields%element_list, &
             particle_in%x(1), particle_in%x(2), & ! inputs
             particle_in%x(1), particle_in%x(2), particle_in%i_elm, &
             particle_in%st(1), particle_in%st(2), ifail) ! outputs
        !p%p = [1.d7,0.]
        !energy = 1.d7 ! 5.12d5 ! Particle energy, including rest energy, in eV
        energy = Etotalini(ini) ! 5.12d5 ! Particle energy, including rest energy, in eV
        !ksi    = 0. !1.d0      ! Cosine of pitch-angle
        ksi    = cosini(ini) !1.d0      ! Cosine of pitch-angle
        !particle_in = p
        particle_out = relativistic_gc_momenta_from_E_cospitch(particle_in,energy,ksi,sim%groups(1)%mass,sim%fields,sim%time)
        particle_in%p = particle_out%p
        write(*,*) 'Initial x(1), x(2), x(3)', particle_in%x
        write(*,*) 'Initial p(1), p(2)', particle_in%p

        ! transform from gc to kinetic:
        chi = dble(ini) ! some gyroangle
        call sim%fields%calc_EBpsiU(sim%time, particle_in%i_elm, &
             particle_in%st, &              
             particle_in%x(3), E, B, psi, U)
        write(6,*) 'HI: B(3):',B
        particle_out_L = relativistic_gc_to_relativistic_kinetic(sim%fields%node_list, &
             sim%fields%element_list, particle_in, sim%groups(1)%mass, B, chi)
        if (particle_out_L%i_elm .eq. 0) cycle

        call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, &
             particle_in%x(1), particle_in%x(2), particle_in%st(1), &
             particle_in%st(2), particle_in%i_elm, particle_out_L%x(1), particle_out_L%x(2), &
             particle_out_L%st(1), particle_out_L%st(2), particle_out_L%i_elm, ifail)

        p = particle_out_L
        p%p = particle_out_L%p
        sim%groups(1)%particles(ini) = particle_out_L
        !sim%groups(1)%particles(ini)%p = particle_out_L%p
        write(6,*) 'HI11: x',sim%groups(1)%particles(ini)%x
        select type (pp=>sim%groups(1)%particles(ini))
              type is (particle_kinetic_relativistic)
              !write(6,*) 'HI11: p',sim%groups(1)%particles(ini)%p
              write(6,*) 'HI11: pp',pp%p
        end select
        write(6,*) 'HI11: particle_out_L%x',particle_out_L%x
        write(6,*) 'HI11: particle_out_L%p',particle_out_L%p
        write(6,*) 'HI11: p%x',p%x
        write(6,*) 'HI11: p%p',p%p
     end select
  enddo ! ini 

!if (i_gc == 0) then
!   do i = 1, n_part
!     chi = dble(i) ! some gyroangle
!     particle_out_L = relativistic_gc_to_relativistic_kinetic(sim%fields%node_list,sim%fields%element_list,particle_in,sim%groups(1)%mass,B,chi)
!     if (particle_out_L%i_elm .eq. 0) cycle
! 
!     call find_RZ_nearby(sim%fields%node_list, sim%fields%element_list, particle_in%x(1), particle_in%x(2), particle_in%st(1), &
!                         particle_in%st(2), particle_in%i_elm, particle_out_L%x(1), particle_out_L%x(2), &
!                         particle_out_L%st(1), particle_out_L%st(2), particle_out_L%i_elm, ifail)
!     
!     !if (ifail .ne. 0) my_ifail = ifail
!   enddo
!endif

else

  write(*,*) '*******************************'
  write(*,*) 'Reading particle restart file'
  call read_simulation_hdf5(sim, 'part_restart.h5')
  write(*,*) 'Time = ', sim%time
  write(*,*) '*******************************' 
  ! Set events to write output data and stop the simulation.
  ! One can use read_jorek_fields_interp_linear or read_jorek_fields_interp_hermite_birkhoff,
  ! and i=-1 (to read jorek_restart.h5 and keep this field at all time) or i=last_file_before_time(sim%time)
  ! (to read a sequel of jorekXXXXX.h5 files and use time-evolving fields)
  events = [event(read_jorek_fields_interp_linear(i=-1)), & 
            event(diag,start=sim%time,step=1d-8),         &
	    event(stop_action(),start=sim%time+1.d-8)]

  ! Run first event to read the JOREK fields
  call with(sim, events, at=0.d0)   

endif

! Check all events conform to the requested timestep
!call check_and_fix_timesteps(timesteps, events)

! Set dpsi/dt=0 (useful e.g. to check the conservation of the total particle energy) 
!call sim%fields%set_flag_dpsidt(.true.)

! Open a file where to write some fields at a given position to test time interpolation routines
!open(22,file='field_vs_t.dat')

! Loop until the simulation is stopped
do while (.not. sim%stop_now)
  ! Extract the next event time
  target_time = next_event_at(sim, events)
  ! Loop over all particle groups
  i_loop: do i=1,1
    ! Compute the number of steps for the particle
    n_steps = nint((target_time - sim%time)/timesteps(i))
	write(*,*) 'n_steps', n_steps
    n_lost = 0

    select type (particles => sim%groups(i)%particles)
    type is (particle_kinetic_relativistic)	
!      !$omp parallel do default(private) &
!      !$omp shared (i, n_steps, timesteps, sim) &
!      !$omp reduction(+:n_lost)	
      j_loop: do j=1,size(particles,1)
        write(6,*) 'particle number ',j 
        k_loop: do k=1,n_steps
          if (particles(j)%i_elm .le. 0) exit

	  sim%time = sim%time + timesteps(i)
	  
! guiding center:   
!          if (i_gc == 1) then  
!           call runge_kutta_fixed_dt_gc_push(sim%fields,sim%time,timesteps(i), &
!               sim%groups(i)%mass,particles(j)) !< push in analytical fields
!             call runge_kutta_fixed_dt_gc_push_jorek(sim%fields,sim%time,timesteps(i), &
!                  sim%groups(i)%mass,particles(j)) !< push in jorek fields
!          endif
! Lorentz:
          if (i_gc == 0) then 
             call volume_preserving_push_jorek(particles(j), sim%fields, sim%groups(i)%mass, &
               sim%time, timesteps(i), ifail)
          endif
!          write(*,*) 'Particle position: ', particles(j)%x(1), particles(j)%x(2), particles(j)%x(3)
!          write(*,*) 'Particle momenta: ', particles(j)%p(1), particles(j)%p(2)

!	  if (modulo(k-1,10000)==0) then	                
!	    call sim%fields%interp_PRZ(sim%time, 1000, [1], 1, 0.5, 0.5, 0.5, P, P_s, P_t, P_phi, P_time, R, R_s, R_t, Z, Z_s, Z_t)	
!	    write(22,'(7e26.16)') sim%time, P, P_time, R, Z
!	  end if

          if (particles(j)%i_elm .le. 0) then
	    n_lost = n_lost + 1
	    write(*,*) 'PARTICLE Nr. ',j,' IS LOST, STOPPING'
	    !stop
            exit k_loop !< or: `cycle j_loop`, should have the same effect
	  end if 		
        end do k_loop !< time steps
      end do j_loop !< particles
!      !$omp end parallel do
    end select
    write(*,*) "number of lost particles: ", n_lost	  
  enddo i_loop !< groups

  ! Update current time and run events
  sim%time = target_time
  call with(sim, events, at=sim%time)
enddo !< event

! Print particle information
write(*,*) 'Final R, Z, phi: ', sim%groups(1)%particles(1)%x(1), sim%groups(1)%particles(1)%x(2), sim%groups(1)%particles(1)%x(3)

! Convert particle to particle_gc to get energy and mu
gyro_angle = 0.
call sim%fields%calc_EBpsiU(sim%time, sim%groups(1)%particles(1)%i_elm, sim%groups(1)%particles(1)%st, &              
                            sim%groups(1)%particles(1)%x(3), E, B, psi, U)
particle_in = sim%groups(1)%particles(1)
call relativistic_gc_to_particle(sim%fields%node_list, sim%fields%element_list, &
                                 particle_in, particle_out_gc, sim%groups(1)%mass, &
				 B, gyro_angle)
write(*,*) 'Final E (eV), mu (eV/T): ', particle_out_gc%E, particle_out_gc%mu

call write_simulation_hdf5(sim, 'part_restart.h5')

! Finalize the simulation
call sim%finalize

end program ex7_jorek

