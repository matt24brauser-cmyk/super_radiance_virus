      program densmat_superradiance
      implicit none
c     constants
      complex, parameter :: zero=cmplx(0.,0.)
      complex, parameter :: eye =cmplx(0.,1.)
      complex, parameter :: one =cmplx(1.,0.)
      real, parameter    :: pi=3.1415926535
c     number of sites
      integer, parameter :: norb  = 10
c     number of spin configurations
      integer, parameter :: nconf = 2**norb
c     Characteristic frequencies
c     radiative time of chromophore [ps]
      real, parameter :: radlife = 200.0 
      real :: gam_rad = 1.0/radlife 
      real :: om_delta = 2.0*pi*10.0 
      real :: om_delta_ini = 0.0
      real :: eta = 1.0/4.0
c     exciton interaction [2pixTHz]
      real :: om_ei = 2.0*pi*1.75
c     exciton interaction-phonon coupling [2pixTHz]
      real :: om_eiph = -2.0*pi*0.1
c     single exciton-phonon coupling strength [2pixTHz]
      real :: om_seph = 2.0*pi*10.0
c     coherent mode frequency [2pixTHz]
      real :: om_coh = 2.0*pi*0.035
c     coherent mode relaxation rate [THz]
      real :: eta_coh = 1.0/65.0 
c     coherent mode inhomogeneous broadening [2pixTHz]
      real :: om_delta_coh = 0.0
c     integration timestep [ps]
      real, parameter    :: dt    = 0.1
c     number of integration timesteps
      integer, parameter :: nstep = 1000
c     number of trajectories
      integer, parameter :: ntraj = 1
c     time-dependent detuning & integral
      real, allocatable    :: domega(:)
c     relaxation exponent
      real, allocatable    :: murelx(:)
c     sigma of noise
      real, allocatable    :: snoise(:)
c     density matrix arrays
      complex, allocatable :: dmat_t(:,:)
c     operator arrays
      complex, allocatable :: si_m(:,:,:)
      complex, allocatable :: s_m(:,:)
      complex, allocatable :: sm_int(:,:) 
      complex, allocatable :: si_d(:,:,:)
      complex, allocatable :: s_d(:,:)
      complex, allocatable :: sd_int(:,:)
      complex, allocatable :: s_pm(:,:)
c     Hamiltonian arrays
      complex, allocatable :: Hev(:,:)
      real, allocatable    :: Hef(:)
      complex, allocatable :: umat(:,:)
c     sampling array
      real, allocatable    :: inten(:)
c     coherent collective mode
      real x_coh, p_coh, pot_coh
      real mu_coh, snoise_coh
c     memory 
      real mem,memmax
c     random numbers
      integer ranseed
      real ran1,rnum,gasdev,mu
c     general
      real maxinten
      integer imax,i,istep,itraj
      real time,hdt
      real smnorm

      write(*,*)'Number of spins',norb
c     number of spin configurations
      write(*,*)'Number of spin configurations',nconf
c     maximum memory 
      memmax=8.0
c     # arrays, 8 bytes per complex, mem in Gb
      mem=((8+2*norb)*8.0*nconf*nconf)/(10.0**9)
      write(*,*)'Memory Gb',mem

      if(mem.lt.memmax)then
        allocate(dmat_t(nconf,nconf),Hev(nconf,nconf),Hef(nconf),
     .           umat(nconf,nconf),sd_int(nconf,nconf),
     .           si_m(nconf,nconf,norb),s_m(nconf,nconf),
     .           si_d(nconf,nconf,norb),s_d(nconf,nconf),
     .           sm_int(nconf,nconf),s_pm(nconf,nconf))
        allocate(domega(norb),murelx(norb),snoise(norb),
     .           inten(nstep))
      else
        write(*,*)'Insufficient memory',mem
        stop
      endif

      write(*,*)'eta          [THz]      = ',eta
      write(*,*)'delta        [2pixTHz]  = ',om_delta
      write(*,*)'Omega e-e    [2pixTHz]  = ',om_ei
      write(*,*)'Omega e-e ph [2pixTHz]  = ',om_eiph
      write(*,*)'Omega e-ph   [2pixTHz]  = ',om_seph
      write(*,*)'Omega coh    [2pixTHz]  = ',om_coh
      write(*,*)'eta coh      [THz]      = ',eta_coh
      write(*,*)'delta coh    [2pixTHz]  = ',om_delta_coh
 
c     initialize the random number generator
      ranseed=1
      rnum=ran1(ranseed)

c     common eta and correlated noise matrix
      do i = 1, norb
        mu = exp(-eta*dt)
        murelx(i) = mu
        snoise(i) = om_delta*sqrt(1.0-mu*mu)
      enddo
      mu_coh = exp(-eta_coh*dt)
      snoise_coh = om_delta_coh*sqrt(1.0-mu_coh*mu_coh)

c     initialize spin operator matrices 
c     si-, S-, |1ix1i|, Sd
      call formSm(norb,nconf,si_m,s_m,si_d,s_d)
c     S+S-
      call cgemm('C','N',nconf,nconf,nconf,one,s_m,nconf,s_m,nconf,
     .           zero,s_pm,nconf)

c     initialize the intensity array
      inten=0.0

c     loop over number of trajectories
      do itraj = 1, ntraj

        if(mod(itraj,10).eq.0)then
          write(*,*)'Trajectory .. ',itraj
        endif
 
c       initialize detuning and integral
c       no initial inhomogeneous broadening due to excitation
        do i = 1, norb
          domega(i) =om_delta_ini*gasdev(ranseed)
        enddo

c       coherent collective mode detuning and integral
        x_coh = 0.0
        p_coh = 0.0

c       initialize density matrix - all spins excited initially
        dmat_t=zero
        dmat_t(nconf,nconf)=one

c       initialize time 
        time=0.0
        hdt=0.5*dt

c       initialize S- interaction representation
        sm_int=s_m
        sd_int=s_d

c       loop over integration timesteps
        do istep = 1, nstep
c        expectation value of S_d at t
         call opexpv(nconf,sd_int,dmat_t,pot_coh) 

c        Euler step - D to t+dt
         call prop_dmat(nconf,(gam_rad*dt),sm_int,dmat_t,smnorm)

c        form H and diagonalize
         call Hdiag(norb,nconf,si_m,s_m,si_d,s_d,s_pm,domega,
     .              (om_seph*x_coh),(om_ei+om_eiph*x_coh),Hef,Hev)
c        form U propagator
         call formU(nconf,dt,Hef,Hev,umat)

c        form S- interaction representation operator at t+dt
         call formSInt(nconf,umat,sm_int)
c        form Sd interaction representation operator at t+dt
         call formSInt(nconf,umat,sd_int)

c        propagate the stochastic detuning to t+dt
         do i = 1, norb
           domega(i)=murelx(i)*domega(i)+snoise(i)*gasdev(ranseed)
         enddo

c        propagate coherent collective mode to t+dt
c        mean-field potential om_seph*<S_d> + om_eiph*(<S+S->-<S_d>)
         pot_coh=om_seph*pot_coh + om_eiph*(smnorm-pot_coh)

c        propagate momentum from t-dt/2 to t+dt/2
         p_coh = p_coh - om_coh*x_coh*dt - pot_coh*dt
         p_coh = mu_coh*p_coh + snoise_coh*gasdev(ranseed)
c        propagate position from t to t+dt
         x_coh = x_coh + om_coh*p_coh*dt

c        sample intensity [unitless - for units of intensity multiply by G_rad*E_rad]
         inten(istep)=inten(istep)+smnorm/real(ntraj)

c        increment time
         time=time+dt

c        print
         if(ntraj.eq.1)then
           write(*,*)time,x_coh,x_coh*0.0025,p_coh,inten(istep)
         endif

        enddo
      enddo

c     intensity maximum and time of maximum
      maxinten=0.0
      imax=1
      do i=1,nstep
        if(inten(i).gt.maxinten)then
          maxinten=inten(i)
          imax=i
        endif
      enddo

c     write results
      open(unit=1,file='maxintensity.txt')
        write(1,*)imax,imax*dt,maxinten
      close(1)
     
      open(unit=1,file='intensity.txt')
      time=0.0
      do istep = 1, nstep
        time=time+dt
        write(1,*) time, inten(istep)
      enddo 
      close(1)

      stop
      end

      subroutine prop_dmat(nconf,dt,s_m,dmat,expv)
      implicit none
      complex, parameter :: zero=cmplx(0.,0.)
      complex, parameter :: eye =cmplx(0.,1.)
      complex, parameter :: one =cmplx(1.,0.)
      complex, parameter :: half=cmplx(0.5,0.)
      integer nconf
      complex s_m(nconf,nconf)
      complex s_pm(nconf,nconf)
      complex dmat(nconf,nconf)
      complex lmat(nconf,nconf)
      complex tmp(nconf,nconf)
      integer i,j
      real dt,expv

c      S-dS+
       call cgemm('N','N',nconf,nconf,nconf,one,s_m,nconf,dmat,nconf,
     .           zero,tmp,nconf)
       call cgemm('N','C',nconf,nconf,nconf,one,tmp,nconf,s_m,nconf,
     .           zero,lmat,nconf)

c      1/2 S+S- 
       call cgemm('C','N',nconf,nconf,nconf,half,s_m,nconf,s_m,nconf,
     .           zero,s_pm,nconf)

c      -( H d + (H d)+ )
       call cgemm('N','N',nconf,nconf,nconf,one,s_pm,nconf,dmat,nconf,
     .           zero,tmp,nconf)

       do i = 1, nconf
         do j = 1, nconf
           lmat(i,j)=lmat(i,j)-tmp(i,j)-conjg(tmp(j,i))
         enddo
       enddo

c      1/2 S+S- d
       call cgemm('N','N',nconf,nconf,nconf,one,s_pm,nconf,dmat,nconf,
     .           zero,tmp,nconf)

c      Tr(S+S- d)
       expv=0.d0
       do i = 1, nconf
         expv=expv+real(tmp(i,i))
       enddo
       expv=2.0*expv

c      Euler step to t+dt
       dmat = dmat + dt * lmat

      return
      end

      subroutine Hdiag(norb,nconf,si_m,s_m,si_d,s_d,s_pm,om,
     .                 om_d,om_int,freq,eigv)
      implicit none
      complex, parameter :: zero=cmplx(0.,0.)
      complex, parameter :: eye =cmplx(0.,1.)
      complex, parameter :: one =cmplx(1.,0.)
      integer norb,nconf
      complex si_m(nconf,nconf,norb)
      complex s_m(nconf,nconf)
      complex s_pm(nconf,nconf)
      complex si_d(nconf,nconf,norb)
      complex s_d(nconf,nconf)
      real om(norb)
      real hamil(nconf,nconf)
      real om_int,om_d
      integer i
      integer info
      real freq(nconf)
      real aux(10*nconf)
      complex eigv(nconf,nconf)

c      form H = om_i si+si- + om_int (S+S--S_d) + om_d S_d
c      om_d = om_seph*x_coh, om_int = om_ei + om_eiph*x_coh
       hamil = 0.0
       do i = 1, norb
         hamil(:,:) = hamil(:,:) + (om(i)+om_d) * real(si_d(:,:,i))
       enddo

c      add interaction term
       hamil(:,:) = hamil(:,:) + om_int * real(s_pm(:,:) - s_d(:,:))

c      diagonalize H
       call ssyev('V','U',nconf,hamil,nconf,freq,aux,10*nconf,info) 

       eigv(:,:)=cmplx(hamil(:,:),0.0)

      return
      end

      subroutine formSInt(nconf,efreq,sm_int)
      implicit none
      complex, parameter :: zero=cmplx(0.,0.)
      complex, parameter :: eye =cmplx(0.,1.)
      complex, parameter :: one =cmplx(1.,0.)
      integer nconf
      complex sm_int(nconf,nconf)
      complex tmp(nconf,nconf)
      complex efreq(nconf,nconf)

c     form U+S-U 
      call cgemm('N','N',nconf,nconf,nconf,one,sm_int,nconf,efreq,nconf,
     .            zero,tmp,nconf) 
      call cgemm('C','N',nconf,nconf,nconf,one,efreq,nconf,tmp,nconf,
     .           zero,sm_int,nconf) 

      return
      end

      subroutine formU(nconf,time,freq,eigv,efreq)
      implicit none
      complex, parameter :: zero=cmplx(0.,0.)
      complex, parameter :: eye =cmplx(0.,1.)
      complex, parameter :: one =cmplx(1.,0.)
      integer nconf
      real time
      complex tmp(nconf,nconf)
      integer i
      real freq(nconf)
      complex efreq(nconf,nconf)
      complex eigv(nconf,nconf)

c      form U = T e^{-ite} T+
       efreq=zero
       do i = 1, nconf
         efreq(i,i) = exp(-eye*time*freq(i))
       enddo
       call cgemm('N','N',nconf,nconf,nconf,one,eigv,nconf,efreq,nconf,
     .            zero,tmp,nconf) 
       call cgemm('N','C',nconf,nconf,nconf,one,tmp,nconf,eigv,nconf,
     .            zero,efreq,nconf)

      return
      end

      subroutine formSm(norb,nconf,si_m,s_m,si_d,s_d)
      implicit none
      complex, parameter :: zero=cmplx(0.,0.)
      complex, parameter :: eye =cmplx(0.,1.)
      complex, parameter :: one =cmplx(1.,0.)
      integer norb,nconf
      complex si_m(nconf,nconf,norb)
      complex s_m(nconf,nconf)
      complex si_d(nconf,nconf,norb)
      complex s_d(nconf,nconf)
      integer isite,iconf,newconf

c     initialize the operators 
      si_m=zero
      s_m=zero
      si_d=zero
      s_d=zero

c     loop over configurations
      do iconf = 0, nconf-1
c        loop over sites
         do isite = 0, norb-1
c          test if excited state
           if (btest(iconf,isite)) then
c            lower the site occupation
             newconf=ibclr(iconf,isite)
c            Si-
             si_m(newconf+1,iconf+1,isite+1)=one
c            Si+Si-
             si_d(iconf+1,iconf+1,isite+1)=one
           endif
         enddo         
      enddo

c     sum Si-
      do isite = 1, norb
        s_m(:,:) = s_m(:,:) + si_m(:,:,isite)
      enddo

c     sum Si+Si-
      do isite = 1, norb
        s_d(:,:) = s_d(:,:) + si_d(:,:,isite)
      enddo

      return
      end

      subroutine opexpv(nconf,op,dmat,expv)
      implicit none
      complex, parameter :: zero=cmplx(0.,0.)
      complex, parameter :: eye =cmplx(0.,1.)
      complex, parameter :: one =cmplx(1.,0.)
      integer nconf
      complex op(nconf,nconf)
      complex dmat(nconf,nconf)
      real expv
      integer iconf
      complex tmp(nconf,nconf)

      call cgemm('N','N',nconf,nconf,nconf,one,op,nconf,dmat,nconf,
     .           zero,tmp,nconf)

c     trace
      expv=0.d0
      do iconf = 1, nconf
        expv=expv+real(tmp(iconf,iconf))
      enddo

      return
      end

      subroutine print_wfn(norb,nconf,psi)
      implicit none
      integer norb,nconf
      complex psi(nconf)
      integer iconf
   
      write(*,*)'Number of sites = ',norb
      do iconf = 0, nconf-1
        write(*,'(i32,B32,2F12.8)')iconf+1,iconf,
     .          real(psi(iconf+1)),aimag(psi(iconf+1))
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c random number generator
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Numerical recipes in Fortran 77 Gaussian distributed random numbers
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
1001    v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1001
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
