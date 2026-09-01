      program virus_particle_hess
      implicit none
c     number of capsomers - pentamers are first Ncaps_pent and then follow from Ncaps_pent+1 hexamers to Ncaps
      integer, parameter :: Ncaps=32
c     number of pentamers 
      integer, parameter :: Ncaps_pent=12
c     number of coat proteins
      integer, parameter :: Ncoat=Ncaps_pent*5+(Ncaps-Ncaps_pent)*6
c     temperature [K]
      real*8, parameter :: temp=300.d0
c     Boltzmann constant [cm-1/K]
      real*8, parameter :: kBoltz=0.695035599492787d0
c     golden ratio
      real*8, parameter :: phi=0.5d0*(1.d0+sqrt(5.d0))
c     virion radius
      real*8, parameter :: rvirion=11.5d0
c     virion icosahedron edge length [nm] (estimated from BMV virion mid-radius of 11.5 nm)
      real*8, parameter :: licosa=rvirion*2.d0/sqrt(phi*phi+1.d0)
c     virion dodecahedron edge length [nm] (hexamers occupy the verteces of a dodecahedron circumscribed in the icosahedron)
      real*8, parameter :: ldodeca=licosa/phi
c     coat protein equilibrium distance [nm] (estimated from capsomers)
      real*8, parameter :: lcoat=2.792d0
c     coat protein mass [Da] (coat protein of BMV) 
      real*8, parameter :: mass_coat_prot=20295.d0
c     conversion factor from Da to system of units defined by energy [cm-1], distance [nm], time [ps]
      real*8, parameter :: mass_conv=1.66053906892d-27*1d+6/1.98645d-23
      real*8, parameter :: m_coatprot=mass_coat_prot*mass_conv
c     nearest neighbor cutoff radius [nm]
      real*8, parameter :: rcut=10.d0
c     binding energy [cm-1] (PNAS paper)
      real*8, parameter :: eps0=12.5d0*kBoltz*300.d0
c     harmonic force constant
      real*8, parameter :: kharm=72.d0*eps0/lcoat**2
c     reduced Planck constant in [cm-1 * ps]
      real*8, parameter :: hbar=5.309188967160051d0
c     Pi
      real*8, parameter :: pi=3.14159265359d0
c     masses
      real*8 mass(Ncoat)
c     positions
      real*8 xyz(3,Ncoat)
      real*8 xyzd(3,Ncoat)
c     Hessian
      integer, parameter :: n3=3*Ncoat
      integer, parameter :: n6=2*n3
      integer, parameter :: lwork=n3*1000
c     classical damping [ps-1] - mass proportional
      real*8, parameter :: eta_a=0.d0
c     classical damping [ps] - frequency proportional
      real*8, parameter :: eta_b=0.8d0
c     Hydrodynamic radius [nm]
      real*8, parameter :: rhydro=0.65d0*lcoat
c     Oseen friction [1/ps] - 0.89 cp [mPa.s=0.001*kg/(m.s)] water viscosity
      real*8, parameter :: oseen_eta=6.d0*pi*0.0019d0*rhydro/0.337d0
      real*8 fricm(n3,n3)
      real*8 fricnm(n3,n3),tmp(n3,n3)
      real*8 hess(n3,n3)
      real*8 hessm(n3,n3)
      real*8 Amat(n6,n6)
      real*8 lre(n6),lim(n6)
      real*8 lre_ord(n6),lim_ord(n6)
      real*8 leig(n6,n6),reig(n6,n6)
      real*8 freq(n3)
      real*8 freq0(n3)
      real*8 aux(lwork)
      real*8 misq(n3)
      integer i,j,info,k,ndeg,imin,ii,jj,kk
      integer imode,istep
      real*8 freqo,eta
      real*8 minlim,ovlp
      real*8 per
      real*8 step
      real*8, parameter :: step_range=1.d+4
      integer, parameter :: nstep=100

      write(*,*)
      write(*,*)'               ------------'
      write(*,*)'             Capsid Hessian program'
      write(*,*)'               ------------'
      write(*,*)

c     initialize positions
      call init_xyz(Ncaps, Ncaps_pent, Ncoat, rvirion, xyz)

c     initialize masses 
      do i = 1, Ncoat
        mass(i)=m_coatprot
        misq(3*(i-1)+1)=1.d0/sqrt(mass(i))
        misq(3*(i-1)+2)=1.d0/sqrt(mass(i))
        misq(3*(i-1)+3)=1.d0/sqrt(mass(i))
      enddo     

c     compute Hessian
      call Hessian(Ncoat,n3,kharm,rcut,xyz,hess)

      do j = 1, n3
        do i = 1, n3
          hess(i,j)=hess(i,j)*misq(i)*misq(j)
        enddo
      enddo

c     diagonalize
      call dsyev('V','U',n3,hess,n3,freq,aux,lwork,info)

      freq0=freq
      do i = 1, n3
        freq(i)=hbar*sqrt(abs(freq(i)))
      enddo

      write(*,*)'#, Degen, L, Freq [cm-1], T [ps]:'
      freqo=freq(7)
      ndeg=0
      k=0
      do i=7,n3
        ndeg=ndeg+1
        if(abs(freq(i)-freqo).gt.0.00001)then
          k=k+1
          eta = eta_a + eta_b * (freqo/hbar)**2
          write(*,'(i4,2x,i3,2x,F3.1,2x,2F14.5)')
     .         k,ndeg-1,real(ndeg-2)/2.0,freqo,2.d0*pi*hbar/freqo 
          ndeg=1
          freqo=freq(i)
        else
          freqo=freq(i)
        endif
      enddo

      write(*,*)
      write(*,*)'hydrodynamic radius [nm] = ',rhydro
      write(*,*)'estimated radius [nm] = ',
     .          sqrt(2.d0*rvirion**2/dble(Ncoat))
      write(*,*)'area ratio = ',rhydro**2/(2.d0*rvirion**2/dble(Ncoat))
      write(*,*)'virion area [nm^2] = ',4.d0*pi*rvirion**2 
      write(*,*)'estimated area [nm^2] = ',dble(Ncoat)*4.d0*pi*rhydro**2
      write(*,*)'friction [ps-1], [ps] = ',oseen_eta,1.d0/oseen_eta
      write(*,*)
 
c     compute the friction matrix
      call Oseen(Ncoat,n3,oseen_eta,rhydro,xyz,fricm)

c     re-compute Hessian
      call Hessian(Ncoat,n3,kharm,rcut,xyz,hessm)

      do j = 1, n3
        do i = 1, n3
          hessm(i,j)=hessm(i,j)*misq(i)*misq(j)
        enddo
      enddo

c     form the A matrix
      Amat=0.d0
      do i = 1, n3
        Amat(i,n3+i)=1.d0
      enddo
      do j = 1, n3
        do i = 1, n3
          Amat(i+n3,j)=-hessm(i,j)
        enddo
      enddo
      do j = 1, n3
        do i = 1, n3
          Amat(i+n3,j+n3)=-fricm(i,j)
        enddo
      enddo

c     diagonalize the A matrix to obtain the T-matrix and l-eigenvalues
      call dgeev('V','V',n6,Amat,n6,lre,lim,leig,n6,reig,n6,
     .           aux,lwork,info)

      write(*,*)'#, Degen, L, Freq [cm-1], T [ps], Damping [ps] ... :'
      freqo=freq(7)
      ndeg=0
      k=0
      do i=7,n3
        ndeg=ndeg+1
        if(abs(freq(i)-freqo).gt.0.00001)then
          k=k+1

          lre_ord=0.d0
          do ii=1,n6
            ovlp=0.d0
            do kk=1,n3
              ovlp=ovlp+hess(kk,i-1)*reig(kk,ii)
            enddo
            lre_ord(ii)=abs(ovlp)
          enddo

          lim_ord=-lre
          do ii = 1, n6-1
            minlim=lre_ord(ii)
            imin=ii
            do jj = ii+1, n6
              if (minlim.ge.lre_ord(jj)) then
                minlim=lre_ord(jj)
                imin=jj
              endif
            enddo
            minlim=lre_ord(ii)
            lre_ord(ii)=lre_ord(imin)
            lre_ord(imin)=minlim
            minlim=lim_ord(ii)
            lim_ord(ii)=lim_ord(imin)
            lim_ord(imin)=minlim
          enddo

          eta = eta_a + eta_b * (freqo/hbar)**2
          write(*,'(i4,2x,i3,2x,F3.1,2x,2F14.5,6F10.3)')
     .         k,ndeg-1,real(ndeg-2)/2.0,freqo,2.d0*pi*hbar/freqo,
     .         (lre_ord(ii),1.d0/lim_ord(ii),ii=n6,n6-2,-1) 

          ndeg=1
          freqo=freq(i)
        else
          freqo=freq(i)
        endif
      enddo

c     inverse mass weight the normal modes       
      do j=1,n3
         do i=1,n3
            hess(i,j)=hess(i,j)*misq(i)
         enddo
      enddo

      imode=161

      open(unit=11,file='mode161.xyz')

      do istep=1, nstep

        step=-step_range*0.5+step_range*(istep-1)/dble(nstep)

c       displace geometry with step along mode imode
        xyzd=xyz
        j=1
        do i = 1, Ncoat
           xyzd(1,i)=xyzd(1,i)+step*hess(j  ,imode)
           xyzd(2,i)=xyzd(2,i)+step*hess(j+1,imode)
           xyzd(3,i)=xyzd(3,i)+step*hess(j+2,imode)
           j=j+3
        enddo

        write(11,*)Ncoat
        write(11,*)
        do i=1, Ncoat
          if (i.le.5*Ncaps_pent) then
            write(11,*)'N  ',0.5*xyzd(1,i),0.5*xyzd(2,i),
     .                       0.5*xyzd(3,i)
          else
            write(11,*)'C  ',0.5*xyzd(1,i),0.5*xyzd(2,i),
     .                       0.5*xyzd(3,i)
          endif
        enddo

      enddo

      close(11)

      return
      end

ccccccccccccccccccccccccccc
c initialize positions
ccccccccccccccccccccccccccc

      subroutine init_xyz(Ncaps, Ncaps_pent, Ncoat, rvirion, xyz_coat)
      implicit none
      integer Ncaps, Ncaps_pent
      real*8 rvirion
      real*8 xyz(3,Ncaps)
      integer i,j,k,kk
      real*8 edgeico,phi,edgedod
      real*8 xyz_coat(3,Ncoat)
      real*8 nvec(3),magn,rvec(3),rmagn,kvec(3)
      real*8 lcoat
      integer Ncoat
      integer hex_type
      real*8, parameter :: pi=3.14159265359d0

        hex_type=2
        xyz=0.d0
        xyz_coat=0.d0

c       golden ratio
        phi=0.5d0*(1.d0+sqrt(5.d0))
c       edge length of the regular icosahedron
        edgeico=rvirion*2.d0/sqrt(phi*phi+1.d0)
c       edge length of the regular dodecahedron
        edgedod=rvirion*2.d0/(phi*sqrt(3.d0))
        lcoat=0.d0

c       initialize the pentamers on the verteces of a regular icosahedron
        xyz(1,1)=phi
        xyz(2,1)=1.d0
        xyz(3,1)=0.d0
        xyz(1,2)=-phi
        xyz(2,2)=1.d0
        xyz(3,2)=0.d0
        xyz(1,3)=phi
        xyz(2,3)=-1.d0
        xyz(3,3)=0.d0
        xyz(1,4)=-phi
        xyz(2,4)=-1.d0
        xyz(3,4)=0.d0

        do i = 1, 4
          xyz(1,i+4)=xyz(3,i)
          xyz(2,i+4)=xyz(1,i)
          xyz(3,i+4)=xyz(2,i)
        enddo

        do i = 1, 4
          xyz(1,i+8)=xyz(2,i)
          xyz(2,i+8)=xyz(3,i)
          xyz(3,i+8)=xyz(1,i)
        enddo

c       scale by the edge length
        xyz(1:3,1:12)=xyz(1:3,1:12)*0.5d0*edgeico

c       initialize the hexamers at the centers of the faces of the icosahedron
c       the centers of faces lie at the verteces of the circumscribed dodecahedron
        xyz(1,13)=phi
        xyz(2,13)=0.d0
        xyz(3,13)=1.d0/phi
        xyz(1,14)=-phi
        xyz(2,14)=0.d0
        xyz(3,14)=1.d0/phi
        xyz(1,15)=phi
        xyz(2,15)=0.d0
        xyz(3,15)=-1.d0/phi
        xyz(1,16)=-phi
        xyz(2,16)=0.d0
        xyz(3,16)=-1.d0/phi

        do i = 13, 16
          xyz(1,i+4)=xyz(2,i)
          xyz(2,i+4)=xyz(3,i)
          xyz(3,i+4)=xyz(1,i)
        enddo

        do i = 13, 16
          xyz(1,i+8)=xyz(3,i)
          xyz(2,i+8)=xyz(1,i)
          xyz(3,i+8)=xyz(2,i)
        enddo

        kk=24
        do i=-1,1,2
          do j=-1,1,2
            do k=-1,1,2
              kk=kk+1
              xyz(1,kk)=dble(i)
              xyz(2,kk)=dble(j)
              xyz(3,kk)=dble(k)
            enddo
          enddo
        enddo

c       scale
        xyz(1:3,13:32)=xyz(1:3,13:32)*0.5d0*phi*edgedod

c       coat protein coordinates
c       pentamers         
        k=0

        if(hex_type.eq.1)then

        do i=1,Ncaps_pent
          do j=Ncaps_pent+1,Ncaps
            nvec(1:3)=xyz(1:3,j)-xyz(1:3,i)
            magn=sqrt(nvec(1)*nvec(1)+nvec(2)*nvec(2)+nvec(3)*nvec(3))
            nvec(1:3)=nvec(1:3)/magn
            if(magn.le.edgedod) then
              k=k+1
              lcoat=magn/3.d0
              xyz_coat(1:3,k)=xyz(1:3,i)+nvec(1:3)*lcoat
            endif
          enddo
        enddo

        else
    
        do i=1,Ncaps_pent
          do j=Ncaps_pent+1,Ncaps
            nvec(1:3)=xyz(1:3,j)-xyz(1:3,i)
            magn=sqrt(nvec(1)*nvec(1)+nvec(2)*nvec(2)+nvec(3)*nvec(3))
            nvec(1:3)=nvec(1:3)/magn
            if(magn.le.edgedod) then
              k=k+1
              lcoat=magn/2.5d0
              rmagn=sqrt(xyz(1,i)**2+xyz(2,i)**2+xyz(3,i)**2)
              rvec(1:3)=xyz(1:3,i)/rmagn
              kvec(1)=rvec(2)*nvec(3)-rvec(3)*nvec(2)
              kvec(2)=rvec(3)*nvec(1)-rvec(1)*nvec(3)
              kvec(3)=rvec(1)*nvec(2)-rvec(2)*nvec(1)
              xyz_coat(1:3,k)=xyz(1:3,i)
     .               +nvec(1:3)*lcoat*0.5d0/tan(pi/5.d0)
     .               +kvec(1:3)*lcoat*0.5d0
            endif
          enddo
        enddo
        
        endif

c       hexamers
        if(hex_type.eq.1)then

        do i=Ncaps_pent+1, Ncaps
          do j=1,Ncaps
            nvec(1:3)=xyz(1:3,j)-xyz(1:3,i)
            magn=sqrt(nvec(1)*nvec(1)+nvec(2)*nvec(2)+nvec(3)*nvec(3))
            nvec(1:3)=nvec(1:3)/magn
            if(magn.le.edgedod*1.1d0.and.i.ne.j) then
              k=k+1
              xyz_coat(1:3,k)=xyz(1:3,i)+nvec(1:3)*lcoat
            endif
          enddo
        enddo

        else

c       hexamers rotated at 30 deg
        do i=1,Ncaps_pent
          do j=Ncaps_pent+1,Ncaps
            nvec(1:3)=xyz(1:3,j)-xyz(1:3,i)
            magn=sqrt(nvec(1)*nvec(1)+nvec(2)*nvec(2)+nvec(3)*nvec(3))
            nvec(1:3)=nvec(1:3)/magn
            if(magn.le.edgedod) then
              k=k+1
              rmagn=sqrt(xyz(1,i)**2+xyz(2,i)**2+xyz(3,i)**2)
              rvec(1:3)=xyz(1:3,i)/rmagn
              kvec(1)=rvec(2)*nvec(3)-rvec(3)*nvec(2)
              kvec(2)=rvec(3)*nvec(1)-rvec(1)*nvec(3)
              kvec(3)=rvec(1)*nvec(2)-rvec(2)*nvec(1)
              xyz_coat(1:3,k)=xyz(1:3,i)
     .               +nvec(1:3)*(magn-lcoat*0.5d0*sqrt(3.d0))
     .               +kvec(1:3)*lcoat*0.5d0
              k=k+1
              xyz_coat(1:3,k)=xyz(1:3,i)
     .               +nvec(1:3)*(magn-lcoat*0.5d0*sqrt(3.d0))
     .               -kvec(1:3)*lcoat*0.5d0
            endif
          enddo
        enddo
        endif

        return

c       print the initial virion coordinates in Ang
        open(unit=1,file='initial_virus.xyz')
        write(1,*)Ncaps
        write(1,*)
        do i=1, Ncaps
          if (i.le.Ncaps_pent) then
            write(1,*)'N  ',xyz(1,i),xyz(2,i),xyz(3,i)
          else
            write(1,*)'C  ',xyz(1,i),xyz(2,i),xyz(3,i)
          endif
        enddo
        close(1)

        open(unit=1,file='initial_virus_coat.xyz')
        write(1,*)Ncoat
        write(1,*)
        do i=1, Ncoat
          if (i.le.5*Ncaps_pent) then
            write(1,*)'N  ',0.5*xyz_coat(1,i),0.5*xyz_coat(2,i),
     .                      0.5*xyz_coat(3,i)
          else
            write(1,*)'C  ',0.5*xyz_coat(1,i),0.5*xyz_coat(2,i),
     .                      0.5*xyz_coat(3,i)
          endif
        enddo
        close(1)

      return
      end

c     elastic network model (ENM) Hessian
      subroutine Hessian(Ncoat,Ncoat3,kharm,rcut,xyz,h)
      implicit none
      integer Ncoat
      integer Ncoat3
      real*8 kharm
      real*8 rcut
      real*8 xyz(3,Ncoat)
      real*8 h(Ncoat3,Ncoat3)
      integer i,j,indi,indj,ii,jj
      real*8 rij2, rcut2
      real*8 xyzij(3)

c     initialize
      h=0.d0

      rcut2=rcut*rcut

c     loop over unique coat-coat protein pairs
      do i = 1, Ncoat

        do j = 1, i-1

c         distance vector
          xyzij(1)=xyz(1,i)-xyz(1,j)
          xyzij(2)=xyz(2,i)-xyz(2,j)
          xyzij(3)=xyz(3,i)-xyz(3,j)
c         distance squared
          rij2=xyzij(1)*xyzij(1)+
     .         xyzij(2)*xyzij(2)+
     .         xyzij(3)*xyzij(3)

c         determine if nearest neighbor
          if (rij2.lt.rcut2) then

           indi=3*(i-1)
           indj=3*(j-1)

           do ii=1,3
             do jj=1,3
               h(indi+ii,indj+jj)=h(indi+ii,indj+jj)
     .                    -kharm*xyzij(ii)*xyzij(jj)/rij2
               h(indj+jj,indi+ii)=h(indj+jj,indi+ii)
     .                    -kharm*xyzij(ii)*xyzij(jj)/rij2
               h(indi+ii,indi+jj)=h(indi+ii,indi+jj)
     .                    +kharm*xyzij(ii)*xyzij(jj)/rij2
               h(indj+ii,indj+jj)=h(indj+ii,indj+jj)
     .                    +kharm*xyzij(ii)*xyzij(jj)/rij2
             enddo
           enddo

          endif

        enddo
      enddo

      return
      end

c     Rotne-Prager friction matrix
      subroutine Oseen(Ncoat,Ncoat3,eta,ra,xyz,h)
      implicit none
      integer Ncoat
      integer Ncoat3
      real*8 eta,ra
      real*8 xyz(3,Ncoat)
      real*8 h(Ncoat3,Ncoat3)
      integer i,j,indi,indj,ii,jj
      real*8 rij2,rij,kscal,dd,dscal,oscal,ra2
      real*8 xyzij(3)

c     initialize
      h=0.d0

      ra2=ra*ra

c     loop over unique coat-coat protein pairs
      do i = 1, Ncoat

        do j = 1, i

          indi=3*(i-1)
          indj=3*(j-1)

          if (j.eq.i) then

c          diagonal term - \delta_ij
           do ii=1,3
              h(indi+ii,indj+ii)=1.d0
           enddo

          else
c          distance vector
           xyzij(1)=xyz(1,i)-xyz(1,j)
           xyzij(2)=xyz(2,i)-xyz(2,j)
           xyzij(3)=xyz(3,i)-xyz(3,j)
c          distance squared
           rij2=xyzij(1)*xyzij(1)+
     .          xyzij(2)*xyzij(2)+
     .          xyzij(3)*xyzij(3)
           rij=sqrt(rij2)

c          non-overlaping spheres
           if(rij.ge.(2.d0*ra))then

           kscal=3.d0*0.25d0*ra/rij
           oscal=2.d0*ra2/rij2
           dscal=kscal*(1.d0+oscal/3.d0)
           oscal=kscal*(1.d0-oscal)

           do ii=1,3
               h(indi+ii,indj+ii)=h(indi+ii,indj+ii)
     .                    +dscal
               h(indj+ii,indi+ii)=h(indj+ii,indi+ii)
     .                    +dscal
           enddo
           do ii=1,3
             do jj=1,3
               h(indi+ii,indj+jj)=h(indi+ii,indj+jj)
     .                    +oscal*xyzij(ii)*xyzij(jj)/rij2
               h(indj+jj,indi+ii)=h(indj+jj,indi+ii)
     .                    +oscal*xyzij(ii)*xyzij(jj)/rij2
             enddo
           enddo

           else
c          overlapping spheres
           dscal=1.d0-9.d0*rij/(32.d0*ra)
           oscal=3.d0/(32.d0*ra)

           do ii=1,3
               h(indi+ii,indj+ii)=h(indi+ii,indj+ii)
     .                    +dscal
               h(indj+ii,indi+ii)=h(indj+ii,indi+ii)
     .                    +dscal
           enddo
           do ii=1,3
             do jj=1,3
               h(indi+ii,indj+jj)=h(indi+ii,indj+jj)
     .                    +oscal*xyzij(ii)*xyzij(jj)/rij
               h(indj+jj,indi+ii)=h(indj+jj,indi+ii)
     .                    +oscal*xyzij(ii)*xyzij(jj)/rij
             enddo
           enddo

           endif

          endif

        enddo
      enddo

      call dmatinv(h,Ncoat3,Ncoat3,dd)

      h=eta*h

      return
      end

      SUBROUTINE DMATINV(A,LDM,N,D)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LDM, N
      REAL*8, INTENT(OUT)   :: D
      REAL*8, INTENT(INOUT) :: A(LDM,*)
      INTEGER             :: I, J, K, L(N), M(N)
      REAL*8                :: BIGA, TEMP
      REAL*8, PARAMETER     :: TOL = 1.0D-12
!
      D = 1.0
!
      DO K = 1, N
        L(K) = K
        M(K) = K
        BIGA = A(K,K)
        DO J = K, N
          DO I = K, N
            IF ( ABS(BIGA).LT.ABS(A(J,I)) ) THEN
              BIGA = A(J,I)
              L(K) = I
              M(K) = J
            END IF
          END DO
        END DO
        J = L(K)
        IF ( J.GT.K ) THEN
          DO I = 1, N
            TEMP = -A(I,K)
            A(I,K) = A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
        I = M(K)
        IF ( I.GT.K ) THEN
          DO J = 1, N
            TEMP = -A(K,J)
            A(K,J) = A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
        IF ( ABS(BIGA).LT.TOL ) THEN
          D = 0.0
          RETURN
        END IF
        DO I = 1, N
          IF ( I.NE.K ) A(K,I) = A(K,I)/(-BIGA)
        END DO
        DO I = 1, N
          DO J = 1, N
            IF ( I.NE.K ) THEN
              IF ( J.NE.K ) A(J,I) = A(K,I)*A(J,K) + A(J,I)
            END IF
          END DO
        END DO
        DO J = 1, N
          IF ( J.NE.K ) A(J,K) = A(J,K)/BIGA
        END DO
        D = MAX(-1.0D25,MIN(1.0D25,D))
        D = D*BIGA
        A(K,K) = 1.0/BIGA
      END DO
!
      K = N
      DO
!
        K = K - 1
        IF ( K.LE.0 ) EXIT
        I = L(K)
        IF ( I.GT.K ) THEN
          DO J = 1, N
            TEMP = A(K,J)
            A(K,J) = -A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
        J = M(K)
        IF ( J.GT.K ) THEN
          DO I = 1, N
            TEMP = A(I,K)
            A(I,K) = -A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
      END DO
!
      END
