! <gde_coagulation.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2024  Matthias Steffen Karl
!*
!*    Contact Information:
!*          Dr. Matthias Karl
!*          Sulzbrackring 13
!*          21037 Hamburg
!*          Germany
!*          email:  mattkar@googlemail.com
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*
!*    The MAFOR code is intended for research and educational purposes. 
!*    Users preparing publications resulting from the usage of MAFOR are 
!*    requested to cite:
!*    1.  Karl, M., Pirjola, L., Grönholm, T., Kurppa, M., Anand, S., 
!*        Zhang, X., Held, A., Sander, R., Dal Maso, M., Topping, D., 
!*        Jiang, S., Kangas, L., and Kukkonen, J., Description and 
!*        evaluation of the community aerosol dynamics model MAFOR v2.0,
!*        Geosci. Model Dev., 15, 
!*        3969-4026, doi:10.5194/gmd-15-3969-2022, 2022.
!*
!*****************************************************************************!
!*    All routines written by Matthias Karl
!* 
!*****************************************************************************!
module gde_coagulation

    use gde_constants,  only            : pi,N_A,k_B
    use gde_constants,  only            : dp
    use gde_constants,  only            : cp_air
    use gde_input_data, only            : NU,AI,AS,CS
    use gde_input_data, only            : MMAX,AMAX  
    use gde_input_data, only            : LAM,VIS
    use gde_init_gas,   only            : rp0,Dfrac
    use gde_sensitiv,   only            : ICOAG
    
    private
   
    public :: coagulation
    public :: coagulation_target
    public :: coagulation_coeff
    public :: collection_kernel

  contains



  subroutine coagulation(temp,DTIME,ROOP,DPA,VPT,N,MASS,IMAX, &
                         IAG,B, COAGS,FLUXM,FLUX)
    !********************************************************************
    !
    !     C  O  A  G  U  L  A  T  I  O  N
    !
    !********************************************************************
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Coagulation module
    !
    !      interface
    !      ---------
    !
    !        input:
    !           IMAX   max. number of bins in mode 
    !           temp   air temperature                [K] 
    !           DTIME  time step                      [s]  
    !           VPT    particle volume in bin         [m^3]
    !           N      particle number conc. in bin   [1/m^3]
    !           ROOP   total particle density         [kg/m^3]
    !           MASS   component mass conc.           [ng/m^3]
    !           DPA    particle diameter in bin       [m]
    !           IAG    coagulation target class, integer matrix
    !           B      coagulation kernel             [m^3/s]
    !
    !        output:
    !           COAGS  coagulation sink               [1/s]
    !           FLUX   collection flux N in bin       [1/m^3/s] 
    !           FLUXM  collection flux m in bin       [1/m^3/s] 
    !
    !      method
    !      ------
    !      Particle from first size bin collides with particles from all 
    !      other size bins. Particle from second size bin collides with 
    !      particles from third to largest size bin. And so on.
    !      An intermediate volume VTOT is calculated
    !
    !      reference
    !      ---------
    !      none
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------
    
     implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX
    REAL( dp), intent(in)                            :: temp,DTIME
    REAL( dp), dimension(MMAX,IMAX),intent(in)       :: VPT,N,ROOP
    REAL( dp), dimension(MMAX,IMAX,AMAX),intent(in)  :: MASS 
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    REAL( dp), dimension(MMAX,MMAX,IMAX,IMAX), intent(in) :: IAG
    REAL( dp), dimension(MMAX,MMAX,IMAX,IMAX),intent(in)  :: B
 
    ! output
    REAL( dp), dimension(MMAX,IMAX),intent(out)      :: FLUX
    REAL( dp), dimension(MMAX,IMAX,AMAX),intent(out) :: FLUXM
    REAL( dp), intent(out)                           :: COAGS
         
    ! local
    REAL( dp) :: VTOTB
    REAL( dp) :: HALF
    REAL( dp) :: FIJK
    REAL( dp) :: FKJK
    integer   :: I,J,L,O,M,MM,A, J1


! Initialize coagulation fluxes
        do M=NU,CS
         do I=1,IMAX
           FLUX(M,I)=0._dp
         end do
        end do

        do M=NU,CS
         do I=1,IMAX
          do A=1,AMAX
           FLUXM(M,I,A)=0._dp
          end do 
         end do
        end do

        COAGS=0._dp


    !------------------------------------------------------------------
    ! The coagulation formulation for the the interacting particles of
    ! size i and j uses the size-splitting operator f_ijk.
    ! Here, f_ijk for the collisional production is:
    !    f_ijk = [v_k+1 - VTOTB) / (v_k+1 - vk)] * HALF 
    !    with HALF = v_k/VTOTB
    !    and VTOTB = v_i + v_j
    !      v_i: VPT(M,O)
    !      v_j: VPT(O,J)
    !      v_k: VPT(O,L)
    ! The size-splitting algorithm gets problematic if volume ratio
    ! becomes > 0.5 and works best with not more than 40 bins
    !------------------------------------------------------------------

! Source terms are updated by the formation of new particles by coagulation.

      do M=NU,CS

          IF (M.EQ.CS) THEN
            MM=CS
          ELSE
            MM=M+1
          ENDIF


          do I=1,IMAX

           do O=M,MM
!MSK 18.11.2024: Modification of J-loop
!        allow production of [O,L] from particles
!        of next lower mode
            !!!do J=I,IMAX
            if (o==m) J1=I
            if (o>m)  J1=1
            do J=J1,IMAX

             !write(6,*) 'B',M,I,B(M,O,I,J)
             VTOTB=VPT(M,I)+VPT(O,J)
             L=INT(IAG(M,O,I,J))


!MSK 01.04.2021 new definition of HALF
!        HALF = v_k/VTOTB = VPT(O,L)/VTOTB
!        except for the first bin where HALF=0
!        and last bin where HALF=1
             IF ((M.EQ.NU).AND.(I.EQ.1)) THEN
                HALF=0.0
             ELSE IF ((M.EQ.CS).AND.(I.EQ.IMAX)) THEN
                HALF=1.0
             ELSE
                HALF= VPT(M,I)/VTOTB
                HALF= HALF*0.5
             ENDIF
!MSK 01.04.2021 end

      !   write(6,'(5I2,6ES12.4)') M,I,O,J,L,VPT(M,I),VPT(O,J),HALF

! BUG-FIX 29.01.2013
! Production Fluxes corrected
!   vpt(o,j) is equal or greater than vpt(m,i)
!   hence production results in new particles in mode o
!   all flux(m,l) were replaced by flux(o,l)
 
             ! number (collisional production)
             if (l == imax) then
               if (o == cs) then
                 flux(o,l)=flux(o,l)+ half*b(m,o,i,j)*n(m,i)*n(o,j)*vtotb/vpt(o,l)
               else
                 FIJK=abs( vpt(o+1,1)-vtotb )
                 FIJK=FIJK/(vpt(o+1,1)-vpt(o,l))
                 IF ((M.EQ.O).AND.(J.EQ.L)) FIJK=1.0  !self-coagulation
                 flux(o,l)=flux(o,l)+half*b(m,o,i,j)*n(m,i)*n(o,j)*FIJK
        !write(6,'(5I3,6ES12.4)') M,I,O,J,L,FIJK,HALF
        if (FIJK.lt.0.0) then
          print*,"1 stop neg f",vpt(o+1,1),vtotb,vpt(o,l),vpt(m,i),vpt(o,j)
          stop
        endif

!MSK 18.11.2024: The condition for FIJK is critical
!        for transfer to the next higher mode
!        changes here must also be done for mass
!MSK 18.11.2024: FIJK not >1.1
!MSK 16.03.2021: FIJK not >2
                 FIJK=vtotb-vpt(o,l-1)
                 FIJK=FIJK/(vpt(o,l)-vpt(o,l-1))
                 !!!FIJK=min(FIJK,2.0)
                 FIJK=min(FIJK,1.1)
                 flux(o+1,1)=flux(o+1,1)+half*b(m,o,i,j)*n(m,i)*n(o,j)*FIJK

        !write(6,'(5I3,6ES12.4)') M,I,O,J,L,FIJK,HALF
        if (FIJK.lt.0.0) then
          print*,"2 stop neg f",vtotb-vpt(o,l-1),(vpt(o,l)-vpt(o,l-1))
          stop
        endif
               endif

             else
               FIJK=vpt(o,l+1)-vtotb
               FIJK=FIJK/(vpt(o,l+1)-vpt(o,l))
               flux(o,l)=flux(o,l)+half*b(m,o,i,j)*n(m,i)*n(o,j)*FIJK
        !write(6,'(5I3,6ES12.4)') M,I,O,J,L,FIJK,HALF
        if (FIJK.lt.0.0) then
          print*,"3 stop neg f",vpt(o,l+1)-vtotb,(vpt(o,l+1)-vpt(o,l))
          stop
        endif
               if (l == 1) then
               FIJK=vtotb-vpt(o-1,imax)
               FIJK=FIJK/(vpt(o,l)-vpt(o-1,imax))
               flux(o,l+1)=flux(o,l+1)+half*b(m,o,i,j)*n(m,i)*n(o,j)*FIJK
         !write(6,'(5I3,6ES12.4)') M,I,O,J,L,FIJK,HALF
         if (FIJK.lt.0.0) then
           print*,"5 stop neg f",vtotb-vpt(o-1,imax),(vpt(o,l)-vpt(o-1,imax))
           stop
        endif

               else
               FIJK=vtotb-vpt(o,l-1)
               FIJK=FIJK/(vpt(o,l)-vpt(o,l-1))
               flux(o,l+1)=flux(o,l+1)+half*b(m,o,i,j)*n(m,i)*n(o,j)*FIJK

         !write(6,'(5I3,6ES12.4)') M,I,O,J,L,FIJK,HALF
         if (FIJK.lt.0.0) then
           print*,"6 stop neg f",vtotb-vpt(o,l-1),(vpt(o,l)-vpt(o,l-1))
           stop
        endif

               endif

             endif

    !    write(6,'(5I3,6ES12.4)') M,I,O,J,L,HALF,b(m,o,i,j),n(m,i),n(o,j)

             !!! mass (collisional production)
             do a=1,amax             
               if (l == imax) then
                 if (o == cs) then
                   fluxm(o,l,a)=fluxm(o,l,a)+ half*b(m,o,i,j)           * &
                        mass(m,i,a)*n(o,j)*vtotb/vpt(o,l)
                 else
                   FIJK=abs( vpt(o+1,1)-vtotb )
                   FIJK=FIJK/(vpt(o+1,1)-vpt(o,l))
                   IF ((M.EQ.O).AND.(J.EQ.L)) FIJK=1.0  !self-coagulation
                   fluxm(o,l,a)=fluxm(o,l,a)+half*b(m,o,i,j)            * &
                        mass(m,i,a)*n(o,j)*FIJK

!MSK 18.11.2024: FIJK not >1.1
!MSK 16.03.2021: FIJK not >2
                   FIJK=vtotb-vpt(o,l-1)
                   FIJK=FIJK/(vpt(o,l)-vpt(o,l-1))
                   !!!FIJK=min(FIJK,2.0)
                   FIJK=min(FIJK,1.1)
                   fluxm(o+1,1,a)=fluxm(o+1,1,a)+half*b(m,o,i,j)        * &
                        mass(m,i,a)*n(o,j)*FIJK

                 endif 
               else
                 FIJK=vpt(o,l+1)-vtotb
                 FIJK=FIJK/(vpt(o,l+1)-vpt(o,l))
                 fluxm(o,l,a)=fluxm(o,l,a)+half*b(m,o,i,j)              * &
                        mass(m,i,a)*n(o,j)*FIJK

                 if (l == 1) then
                 FIJK=vtotb-vpt(o-1,imax)
                 FIJK=FIJK/(vpt(o,l)-vpt(o-1,imax))
                 fluxm(o,l+1,a)=fluxm(o,l+1,a)+half*b(m,o,i,j)          * &
                        mass(m,i,a)*n(o,j)*FIJK
                 else
                 FIJK=(vtotb-vpt(o,l-1))
                 FIJK=FIJK/(vpt(o,l)-vpt(o,l-1))
                 fluxm(o,l+1,a)=fluxm(o,l+1,a)+half*b(m,o,i,j)          * &
                        mass(m,i,a)*n(o,j)*FIJK

                 endif
               endif
             ! END BUG-FIX 29.01.2013               
             end do
            end do
           end do
          end do


! Source terms are updated by removal of particles by coagulation.
    !------------------------------------------------------------------
    ! The size-splitting operator for collisional loss is (1 - f_kjk)
    ! Here, f_kjk for the collisional loss is:
    !    f_kjk = [v_k+1 - VTOTB) / (v_k+1 - vk)] * HALF 
    !    with HALF = v_k/VTOTB
    !    and VTOTB = v_k + v_j
    !      v_j: VPT(O,J)
    !      v_k: VPT(M,I)
    !------------------------------------------------------------------
          do I=1,IMAX
           do O=M,MM
            do J=1,IMAX

              VTOTB=VPT(M,I)+VPT(O,J)
              IF ((M.EQ.CS).AND.(I.EQ.IMAX)) THEN 
                FKJK = 1.0
              ELSE
                FKJK = VPT(M,I)/VTOTB
              ENDIF
              HALF =1.0 - FKJK

       !   write(6,'(4I2,6ES12.4)') M,I,O,J,VPT(M,I),VPT(O,J),HALF

              ! number (collisional loss)

              FLUX(M,I)=FLUX(M,I)- HALF* B(M,O,I,J)*N(M,I)*N(O,J)


             !! Constraint: coagulation flux not greater than number conc. allows
             !! 28.04.2018 needs negative sign
              IF ((abs(FLUX(M,I)*DTIME) .GT. N(M,I)) .AND. (FLUX(M,I) .LT. 0._dp)) THEN
                !write(6,*) 'warning',M,O,I,J,FLUX(M,I),N(M,I),B(M,O,I,J),N(O,J)
                 FLUX(M,I)= (-1.)*N(M,I)/DTIME 
                write(6,*) 'loss flux gt N'
                !stop
              ENDIF

              !!! mass (collisional loss)
              do A=1,AMAX

                 FLUXM(M,I,A)=FLUXM(M,I,A) - HALF*  B(M,O,I,J)*MASS(M,I,A)*N(O,J)


               IF ((abs(FLUXM(M,I,A)*DTIME).GT.MASS(M,I,A)).AND.(FLUXM(M,I,A).LT.0.)) THEN
                 FLUXM(M,I,A)=MASS(M,I,A)/DTIME
               ENDIF
              end do               
            end do
          end do
           ! coagulation sink 
            COAGS=COAGS+(B(1,M,1,I)*N(M,I))          
          end do
      end do
  
  
  end subroutine coagulation


  subroutine coagulation_target(VPT,IMAX,IAG)
    !----------------------------------------------------------------------
    !
    !****  Compute initial coagulation target classes
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      compute initial coagulation target classes IAG
    !
    !      interface
    !      ---------
    !
    !        input:
    !           IMAX   max. number of bins in mode 
    !           VPT    particle volume in bin         [m^3]
    !
    !        output:
    !           IAG    coagulation target class, integer matrix
    !
    !      method
    !      ------
    !      compute initial coagulation target classes IAG
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      none
    !
    !-----------------------------------------------------------------
     
     implicit none   
    
    ! input
    INTEGER, intent(in)                                   :: IMAX
    REAL( dp), dimension(MMAX,IMAX),intent(in)            :: VPT
    ! output
    REAL( dp), dimension(MMAX,MMAX,IMAX,IMAX), intent(out):: IAG

    REAL( dp) :: VTOTB   
    integer   :: I,J,K,L,O,M,MM
    
      L=1 !dummy initialization 
      do M=NU,CS
       MM=M
    ! Compute coagulation target classes
       do I=1,IMAX
        do O=MM,CS
         do J=1,IMAX
           VTOTB=VPT(M,I)+VPT(O,J)
           IF (VTOTB .GE. VPT(O,IMAX)) THEN
             L=IMAX
           ELSE
            do K=J,IMAX-1
              IF (VTOTB .GE. VPT(O,K) .AND. VTOTB .LT. VPT(O,K+1)) THEN
                L=K
                exit
              ENDIF
             end do
           ENDIF
           IAG(M,O,I,J)=L
         end do
        end do
       end do
      end do
      
    
    end subroutine coagulation_target 

  subroutine coagulation_coeff(temp,ROOP,DPA,coag,IMAX)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      calculation of coagulation coefficient B
    !
    !      interface
    !      ---------
    !           IMAX   max. number of bins in mode 
    !           temp   air temperature                [K] 
    !           ROOP   total particle density         [kg/m^3]
    !           DPA    particle diameter in bin       [m]
    !
    !        output:
    !           coag   coagulation kernel             [m^3 s^-1]
    !
    !      method
    !      ------
    !      calculates coagulation coefficients
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      none
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    INTEGER, intent(in)                              :: IMAX
    REAL( dp), intent(in)                            :: temp
    REAL( dp), dimension(MMAX,IMAX),intent(in)       :: ROOP
    REAL( dp), dimension(MMAX,0:(IMAX+1)),intent(in) :: DPA
    REAL( dp), dimension(MMAX,MMAX,IMAX,IMAX), intent(out):: coag

    real( dp)  :: AA, KNA, CC
    real( dp)  :: CIJ,GIJ,DIJ,DIFIJ
    real( dp)  :: rp0m
    real( dp)  :: KNP,EPS,nspher,rarea,ra1,rmobil
    real( dp)  :: AH
    real( dp)  :: ru, rui,pr,sr,dr
    real( dp)  :: a,b
    real( dp)  :: z1,z2,wc,wk,ve

    real( dp), dimension(MMAX,IMAX)           :: C,DIF,G,LL
    real( dp), dimension(MMAX,0:(IMAX+1))     :: DPAC
    real( dp), dimension(MMAX,0:(IMAX+1))     :: DPAM

    real( dp), dimension(2,40)                :: ir
    real( dp), dimension(40)                  :: r

    ! Lemmetty et al. (2008)
    !real( dp), parameter :: sootfrac = 1.00     ! soot mass fraction
    !real( dp), parameter :: rp0      = 2.50e-9  ! radius of primary spherules [m]
    !real( dp), parameter :: Dfrac    = 2.50     ! fractal dimension           [-]
    ! Jacobson and Seinfeld (2004)
    real( dp), parameter :: sootfrac = 0.50      ! soot mass fraction
    !real( dp), parameter :: rp0      = 13.50e-9  ! radius of primary spherules [m]
    !real( dp), parameter :: Dfrac    = 1.70      ! fractal dimension           [-]
    
    integer :: I,J,M,O
    integer :: n,ii

    ! convert radius of primary spherules from [nm] to [m]
    if ( (ICOAG.eq.2).or.(ICOAG.eq.5) ) rp0m=rp0*1.e-9

    ! first calculate C(M,I) and G(M,I) for all bins
    do M=1,MMAX
     do I=1,IMAX

       if ( (ICOAG.eq.2).or.(ICOAG.eq.5) ) then
       ! Collision radius becomes fractal radius rc rf
          nspher = ( 4*pi*(DPA(M,I)*0.5)**3./3 ) /            &
                   ( 4*pi*rp0m**3./3 )
          DPAC(M,I) = 2.* rp0m * nspher**(1./Dfrac)
          DPAC(M,I) = (1.0-sootfrac)*DPA(M,I)                 &
                    + sootfrac* DPAC(M,I)
          !write(6,*) 'dp_corr',m,i,DPA(M,I),DPAC(M,I)
       ! Mobility radius rm for Knudsen number and diffusion coeff
          ra1   = min( (1.0+0.67_dp*(nspher-1.0_dp) ),   &
                       (Dfrac*nspher**(2.0_dp/Dfrac)/3) )
          rarea = rp0m * sqrt( max( (nspher**(2.0_dp/3)),ra1 ) )
          rmobil = max( (DPAC(M,I)*0.5/(log(DPAC(M,I)*0.5/rp0m)+1.0_dp) ), &
                        (DPAC(M,I)*0.5*((Dfrac-1.0_dp)/2.0_dp)**0.7_dp) )
          DPAM(M,I) = 2.0_dp*max( rmobil,rarea )
          if ( DPA(M,I)*0.5_dp .lt. rp0m ) then
            DPAM(M,I) = DPAC(M,I)
          endif
        !write(6,*) 'dp_mobi',m,i,rp0m,Dfrac,DPA(M,I),DPAC(M,I),DPAM(M,I)
       else
          DPAC(M,I) = DPA(M,I)
          DPAM(M,I) = DPA(M,I)
       endif

       !!! MP(M,I) from dens*vol  (does not work:)
       AA=(1./6.)*ROOP(M,I)*pi*(DPAM(M,I))**3.  ! AA in [kg] 
       !!! THIS IS AN IMPORTANT CONSTRAINT   
       !AA=max(MP(M,I),1.e-22) 
       !!! DO NOT CHANGE
       
       ! Thermal velocity C [m s^-1]
       ! Boltzmann constant here in [J/K] == [kg m^2 s^-2 K^-1]
       ! particle mass in [kg]
       C(M,I)=SQRT(8._dp*k_B*temp/(pi*AA))
       !write(6,*) M,I,C(M,I),AA
       ! Cunnningham slip-flow correction [-]
       ! CC = 1 + Kn_a(A + B*exp(-C/Kn_a))
       ! with A,B,C from Rogak and Flagan [J. Coll. Interface Sci., 151, 1992]
       ! Kn_a is Knudsen number of air 
       KNA = 2._dp*LAM/DPAM(M,I)
       CC=1._dp+KNA*(1.257_dp+0.4_dp*exp(-1.1_dp/KNA))
       ! Particle Diffusion Coefficient   [m^2 s^-1]
       DIF(M,I)=(k_B*temp*CC/(3.*pi*VIS*DPAM(M,I)))
       ! Particle mean free path [m]
       LL(M,I)=8._dp*DIF(M,I)/(pi*C(M,I))
       ! Mean distance from center of a sphere [m]
       G(M,I) = (  ( (DPAM(M,I)+LL(M,I))**3.                   &
               - (DPAM(M,I)**2.+LL(M,I)**2.)**1.5 )            &
               / (3.*DPAM(M,I)*LL(M,I))  )                     &
               - DPAM(M,I)                
       !write(6,*) M,I,LL(M,I),DIF(M,I),G(M,I)
       !write(6,*) M,I,ROOP(M,I),AA,C(M,I),CC
     end do
    end do


    ! second calculate coag(M,O,I,J) for all bins [m^3 s^-1]
    do M=1,MMAX
     do I=1,IMAX
      do O=1,MMAX
       do J=1,IMAX
         DIJ=DPAC(M,I)+DPAC(O,J)
         DIFIJ=DIF(M,I)+DIF(O,J)
         CIJ=SQRT(C(M,I)**2._dp+C(O,J)**2._dp)
         GIJ=SQRT(G(M,I)**2._dp+G(O,J)**2._dp)
         ! ICOAG.eq.3 van der Waals and viscous forces correction
         if (ICOAG.eq.3) then
           !Particle pair Knudsen number
           KNP=SQRT((LL(M,I)**2._dp + LL(O,J)**2._dp))      &
                  /(0.5*DPAC(M,I)+0.5*DPAC(O,J))
           if (KNP.lt.0.1)  EPS=1.0
           if (KNP.ge.0.1 .and. KNP.le.10.0)  EPS=0.4_dp*LOG(KNP)+2.0_dp
           if (KNP.gt.10.0)  EPS=3.0
           !write(6,*) 'knp eps',m,i,o,j,KNP,EPS
         else
           EPS=1.0
         endif
         coag(M,O,I,J)=2.*pi*DIJ*DIFIJ *EPS                 &
                  /(DIJ/(DIJ+2._dp*GIJ)+8.*DIFIJ/(CIJ*DIJ))
         !write(6,*) 'incoag',M,O,I,J,coag(M,O,I,J),C(M,I),C(O,J)  
         !write(6,*) 'incoag',M,I,O,J,coag(M,O,I,J)


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! ICOAG.eq.4 explicit correction for van-der-Waals/viscous forces
         ! M.Z. Jacobson (2005), Fundamentals of Atmospheric Modelling,
         !  Second Edition, page 513
         ! M.Z. Jacobson & J.H. Seinfeld (2004), Evolution of nanoparticle
         !  size and mixing state near the point of emission,
         !  Atmos. Environ. 38, 1839-1850, doi:10.1016/j.atmosenv.2004.01.014
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! ve is the van-der-Waals/viscous collision correction factor
         ! The integrals are approximated by numerical integration using
         ! Gauss-Legendre Quadrature.
         ! the integral is approximated by:
         ! int(a,b) f(x)dx = (b-a)/2 * sum(i,n) wi*f( (a+b)/2 + xi*(b-a)/2 )
         if ( (ICOAG.eq.4).or.(ICOAG.eq.5) ) then
           ru=DPAC(M,I)*0.5_dp
           rui=DPAC(O,J)*0.5_dp
           pr=ru*rui
           sr=ru+rui
           dr=ru-rui

           !Hamaker constant of water AH/(kB*T)=20
           AH=20._dp*k_B*temp
           !print *,"Hamaker AH:",AH  !(7.89698e-20)
           !number of gaussian evaluation points = 20
           !n = 40
           !MSK 2021-03-07 speeed up the evaluation
           n = 20
           !limits of integration
           a = 0._dp
           b = 1._dp/(1._dp+rui/ru)
           ir = gaussquad(n)         
           !function to integrate (here simple exp(x) )
           !ir(1,:) is xi and ir(2,:) is wi
           !z = (b-a)/2*dot_product(ir(2,:),exp((a+b)/2+ir(1,:)*(b-a)/2))

           do ii = 1,n
             r(ii) = (a+b)/2 + ir(1,ii)*(b-a)/2
           enddo

           !wc: Correction factor for the continuum regime
           !Van-der-Waals interaction potential Ep and
           !viscous force correction factor Dratio
           !integral:
           !z1=int (Dratio(r)*exp(Ep(r))*1/r^2) dr

           z1= (b-a)/2*dot_product( ir(2,:),            &
               !Dratio(r)
               ( 1._dp + (2.6_dp*pr/(sr)**2)*                 &
               !value in sqrt must not be negative
               !    sqrt(pr/(sr*(r(:)-ru-rui)))  + &
                 sqrt( abs(pr/(sr*(r(:)-ru-rui))) )    +      & 
                 pr/(sr*(r(:)-ru-rui)) )*                     &
               !exp(Ep)
               exp( (-AH/6*k_B*temp)* ( 2*pr/(r(:)**2-sr**2)+ &
                  2*pr/(r(:)**2-dr**2)  +                     &
                  log((r(:)**2-sr**2)/(r(:)**2-dr**2)) ) )*   &
               !1/r^2
               1/r(:)**2           )
           
           !wk: Correction factor for the free-molecular regime
           !integral:
           !z2=int ( (dEp(r)/dr+r*d^2Ep(r)/dr^2)*  &
           !        exp(-((r/2)*dEp(r)+Ep(r))*r^2 ) dr
           !Differentials dEp(r)/dr and d^2Ep(r)/dr^2 obtained
           !with Matlab symbolic package
           
           z2= (b-a)/2*dot_product( ir(2,:),            &
               !dEp(r)/dr
                 (  AH*r(:)*( 2*pr*(r(:)**2-sr**2)**2+        &
                 2*pr*(r(:)**2-dr**2)**2             +        &
                 (r(:)**2-sr**2)*(r(:)**2-dr**2)     *        &
                 (dr**2-sr**2) )                     /        &
                 ( 3._dp*(r(:)**2-sr**2)**2 *                 &
                 (r(:)**2-dr**2)**2 )  +                      &    
               !r*d^2Ep(r)/dr^2
                 r(:)*AH *( -8._dp*pr*r(:)**2          *      &
                 (r(:)**2-sr**2)**3                    -      & 
                 8._dp*pr*r(:)**2*(r(:)**2-dr**2)**3   +      &
                 2*pr*(r(:)**2-sr**2)**3               *      &
                 (r(:)**2-dr**2)                       +      &
                 2*pr*(r(:)**2-sr**2)*(r(:)**2-dr**2)**3 +    &
                 2*r(:)**2*(r(:)**2-sr**2)**2          *      &
                 (r(:)**2-dr**2)*(dr**2-sr**2)         +      &
                 2*r(:)**2*(r(:)**2-sr**2)             *      &
                 (r(:)**2-dr**2)**2                    *      &
                 (sr**2-dr**2)                   )& !      +      &
                 !from here mixed terms (cause oscillation)
                 !(r(:)**2-sr**2)**2*(r(:)**2-dr**2)    *      &
                 !( -4._dp*r(:)**2*(r(:)**2-sr**2)      -      &
                 !  (r(:)**2-dr**2)**2 + (r(:)**2-dr**2) *     &
                 !  (5._dp*r(:)**2-sr**2) )   )                &       
                   /       &
                 ( 3._dp*(r(:)**2-sr**2)**3 *                 &
                 (r(:)**2-dr**2)**3 )                         &
                 )  *                                         &
               !exp(-((r/2)*dEp(r)+Ep(r))
                 exp(   (-1._dp/k_B*temp)*(  (r(:)/2) *       &
                 ( AH*r(:)*( 2*pr*(r(:)**2-sr**2)**2 +        &
                 2*pr*(r(:)**2-dr**2)**2             +        &
                 (r(:)**2-sr**2)*(r(:)**2-dr**2)     *        &
                 (dr**2-sr**2) )                     /        &
                 ( 3._dp*(r(:)**2-sr**2)**2 *                 &
                 (r(:)**2-dr**2)**2 ) ) )  +                  &
                 (-AH/6._dp)*( 2*pr/(r(:)**2-sr**2) +         &
                  2*pr/(r(:)**2-dr**2)  +                     &
                  log((r(:)**2-sr**2)/(r(:)**2-dr**2))  )     &
                 )*                                           &
               !r^2
                 r(:)**2                  )

           !calculate the collision correction factor ve
           wc=1._dp/(sr*z1)
           wk=(-1._dp/(2*sr**2*k_B*temp))*z2
           !scale by 1/20 for more realistic value
           wk=wk*0.05_dp
           ve=wc *(1._dp+(4._dp*DIFIJ/(CIJ*(ru+rui))) )
           ve=ve/( 1._dp+(wc/wk)*(4._dp*DIFIJ/(CIJ*(ru+rui))) )

           ve=max(ve,1.0_dp)

           !modification of Brownian kernel by van-der-Waals/viscous
           coag(M,O,I,J)=coag(M,O,I,J)*ve

           !Test with Particle pair Knudsen number
           !KNP=SQRT((LL(M,I)**2._dp + LL(O,J)**2._dp))      &
           !       /(0.5*DPAC(M,I)+0.5*DPAC(O,J))
           !if ( (M==2).and.(I==10) ) then
           !  print *,'KNP',M,I,O,J,KNP,wc,wk,ve
           !endif
         endif

        ! if ( (DPA(M,I).gt.9.6e-09) .and. (DPA(M,I).lt.10.2e-09) ) then
        !   write(6,*) 'beta',M,I,O,J,DPA(O,J),DPAC(O,J),coag(M,O,I,J)
        !   !write(6,*) 've',M,I,O,J,wk,wc,ve
        ! endif

       end do
      end do     
     end do
    end do


  end subroutine coagulation_coeff


  subroutine collection_kernel(IMAX,temp,press,DPA,vterm,    &
                               DPcrit,KCOA,  KCOL )
    !----------------------------------------------------------------------
    !
    !****  COLLECTION KERNEL
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      The collection kernel K describes the interaction
    !      of two colliding particles. The collision partners
    !      are assummed to fall with their terminal velocities
    !      Vinf. The smaller collision partner (s) will be 
    !      collected by the larger collision partner (b) as soon
    !      as it is inside the swept volume of the collector
    !      KCOL = Ecoa*Ecol*|Vinf_B - Vinf_s|*pi*(RB+rs)**2
    !      Ecoa: coalescence efficieny, Table 4, Kirkweg (2003)
    !      Ecol: collision efficiency, Table 2, Kirkweg (2003)
    !
    !
    !      interface
    !      ---------
    !
    !        input:
    !           temp     air temperature                [K]
    !           press    air pressure                   [Pa]
    !           DPA      wet particle diameter          [m]
    !           VPT      particle volume                [m^3]
    !           ROOPW    particle density               [kg/m^3]
    !           vterm    terminal velocity              [m/s]
    !           DPcrit   critical particle diameter     [m]
    !           KCOA     coagulation coefficient        [m^3/s]
    !
    !
    !        output:
    !           KCOL     collection kernel              [m^3/s]
    !
    !
    !      method
    !      ------
    !      Collision/Coalescence as function of size
    !      of the collision partners
    !      For mode AI to CS (excluding NU mode)
    !
    !      reference
    !      ---------
    !      Ventilation coefficient:
    !      J.M. Straka,
    !      Cloud and Precipitation Microphysics,
    !      Principles and Paramterizations,
    !      Cambridge University Press, Cambridge, UK,
    !      p.116-118, 2009.
    !
    !      Collision efficiency:
    !      K. Young,
    !      Microphysical Cloud Processes,
    !      Oxford University Press, New York, USA, 1993.
    !
    !      Coalescence efficiency:
    !      K.V. Beard and H.T. Ochs
    !      Collection and Coalescence Efficiencies 
    !      for Accretion,
    !      J. Geophys. Res., 89, D5, 7165-7169, 1984.
    !
    !      Collection kernel:
    !      A. Kerkweg, S. Wurzler, T. Reisin, A. Bott,
    !      On the cloud processing of aerosol particles:
    !      An entraining air-parcel model with two-dimensional
    !      spectral cloud microphysics and a new formulation
    !      of the collection kernel,
    !      Q.J.R. Meteorolo. Soc., 129, 1-19, 2003
    !
    !      modifications
    !      -------------
    !      none
    !
    !------------------------------------------------------------------

    implicit none

    ! input
    integer, intent(in)                                :: IMAX
    real( dp), intent(in)                              :: temp      ! [K]
    real( dp), intent(in)                              :: press     ! [Pa]
    real( dp), dimension(MMAX,0:(IMAX+1)), intent(in)  :: DPA       ! [m]
    real( dp), dimension(MMAX,IMAX)                    :: vterm     ! [m/s]
    real( dp), intent(in)                              :: DPcrit    ! [m]
    real( dp), dimension(MMAX,MMAX,IMAX,IMAX), intent(in)  :: KCOA  ! [m^3/s]

    ! output
    real( dp), dimension(MMAX,MMAX,IMAX,IMAX), intent(out) :: KCOL  ! [m^3/s]


    ! local
    ! specific gas constant for dry air [J/(kgK)]
    real( dp), parameter                               :: R=287.058
    real( dp),dimension(MMAX,MMAX,IMAX,IMAX)           :: ECOL
    real( dp),dimension(MMAX,MMAX,IMAX,IMAX)           :: ECOA
    real( dp), dimension(MMAX,IMAX)                    :: DIF
    real( dp), dimension(MMAX,IMAX)                    :: VEN

    real( dp)                                          :: rho_air
    real( dp)                                          :: MYY
    real( dp)                                          :: NYY
    real( dp)                                          :: KNA
    real( dp)                                          :: CC
    real( dp)                                          :: DIFFCO
    real( dp)                                          :: NRE
    real( dp)                                          :: NSC

    real( dp)                                          :: rc
    real( dp)                                          :: rlim
    real( dp)                                          :: rs
    real( dp)                                          :: RB
    real( dp)                                          :: difs
    real( dp)                                          :: vterms
    real( dp)                                          :: acoa
    real( dp)                                          :: bcoa
    real( dp)                                          :: xcoa
    real( dp)                                          :: betacoa

    integer                                            :: I,M,J,O
    integer                                            :: MM,K


      ! Initialisation
      KCOL(:,:,:,:)   = 0._dp
      ECOL(:,:,:,:)   = 0._dp
      ECOA(:,:,:,:)   = 0._dp
      DIF(:,:)        = 0._dp
      VEN(:,:)        = 1._dp


      ! Critical radius
      rc      = DPcrit*0.5_dp
      rlim    = 4.00e-7_dp      ! 400 nm
      rlim    = MIN(rlim,rc)

      ! Air density
      ! [kg/m^3] 1.2928
      rho_air = press / R / temp      

      ! Dynamic viscosity of air
      ! [kg/m/s]
      MYY= (1.832e-5_dp*(temp**(1.5_dp))*406.4_dp)/            &
          (5093._dp*(temp+110.4_dp))

      ! Kinematic viscosity of air
      ! [m^2/s]
      NYY     = MYY / rho_air


      ! First calculate in all bins
      ! Particle diffusion coefficient
      ! for Brownian Diffusion DIF
      ! and Ventilation coefficient VEN

      do M=AI,CS
        do I=1,IMAX

          ! Cunnningham slip-flow correction [-]
          ! CC = 1 + Kn_a(A + B*exp(-C/Kn_a))
          ! with A,B,C from Rogak and Flagan 
          ! [J. Coll. Interface Sci., 151, 1992]
          ! LAM and VIS in gde_input_data
          ! Kn_a is Knudsen number of air 
          KNA = 2._dp*LAM / DPA(M,I)
          CC  = 1._dp+KNA*(1.257_dp+0.4_dp*exp(-1.1_dp/KNA))

          ! Particle Diffusion Coefficient   [m^2/s]
          DIF(M,I) = ( k_B * temp * CC / (3.*pi*VIS*DPA(M,I)) )

          ! Vapor diffusivity [m^2/s]
          DIFFCO  = (k_B*temp*CC) /                            &
                    (6._dp*pi*MYY*0.5_dp*DPA(M,I))

          ! Reynold number [-]
          NRE = (vterm(M,I)*DPA(M,I)) / NYY

          ! Schmidt number [-]
          NSC = NYY / DIFFCO

          ! Ventilation coefficent [-]
          ! Two conditions depending on NRE and NSC
          if (NSC**(1./3.)*sqrt(NRE) .lt. 1.4) then
            VEN(M,I) = 1._dp + 0.108_dp*                       &
                       (NSC**(1./3.)*sqrt(NRE))**2.
          else
            VEN(M,I) = 0.78_dp + 0.308_dp*                     &
                       NSC**(1./3.)*sqrt(NRE)
          endif

      ! print *,'DpA CASE DIF',M,I,DPA(M,I),NSC**(1./3.)*sqrt(NRE),VEN(M,I),DIF(M,I)

         enddo       ! I
      enddo         ! M


      ! Second calculate the collection kernel KCOL
      ! KCOL = Ecoa*Ecol*|Vinf_B - Vinf_s|*pi*(RB+rs)**2
      !  Ecoa: coalesence efficieny, Table 4, Kirkweg (2003)
      !  Ecol: collision efficiency, Table 2, Kirkweg (2003)
      ! Accretion is the growth of large partner by 
      ! sweeping the smaller partner.
      ! Large collision partner X(O,J)
      !  with radius RB
      ! Small collision partner X(M,I)
      !  with radius rs
      ! Outer loop: small partner
      ! Inner loop: large partner
      ! Coalescence:
      !   small partner rs: 1.0 - 30 um
      !   large partner RB: >50 um
      !   for small partner rs < 1um -> Ecoa=1
      !   for large partner RB < 50 um -> Ecoa=1

      do M=NU,CS
        do I=1,IMAX

          rs     = DPA(M,I)*0.5_dp
          difs   = DIF(M,I)
          vterms = vterm(M,I)

          do O=AI,CS
            do J=1,IMAX

              RB      = DPA(O,J)*0.5_dp
              betacoa = 0._dp
              bcoa    = 0._dp
              acoa    = 0._dp
              xcoa    = 0._dp

              IF (RB.GT.rs) THEN

                ! Collision efficiency [-]
                ! Ecol = 4*RB*DIF*VEN /
                !        (rs+RB)**2*|Vinf_B-Vinf_s|

                ECOL(M,O,I,J) = 4._dp*RB*difs*VEN(O,J)       /  &
                      ( (RB+rs)**2.*ABS(vterm(O,J)-vterms) )

                ! Coalescence efficiency [-]

                ! RB >= 50 um
                if (RB.ge.5.e-5) then
                  if ((rs.ge.1.e-6).and.(rs.le.3.e-5)) then

                    betacoa = LOG(rs*1.e6)                   +  &
                          0.44_dp*LOG(RB*1.e6/200._dp)
                    bcoa    = 0.0946_dp * betacoa - 0.319_dp
                    acoa    = SQRT(bcoa**2.0 + 0.00441_dp)
                    xcoa    = (acoa-bcoa)**(1._dp/3._dp)     -  &
                            (acoa+bcoa)**(1._dp/3._dp)
                    ECOA(M,O,I,J) = xcoa + 0.459_dp

                  else if (rs.gt.3.e-5) then  ! rs > 30 um
                    ECOA(M,O,I,J) = 0.6_dp
                  else                        ! rs < 1 um 
                    ECOA(M,O,I,J) = 1._dp
                  endif

                ! RB < 50 um
                else
                  if (rs.gt.3.e-5) then       ! rs > 30 um
                    ECOA(M,O,I,J) = 0.6_dp
                  else 
                    ECOA(M,O,I,J) = 1._dp
                  endif
                endif

              ! Collection kernel [m^3/s]

                KCOL(M,O,I,J) = ECOA(M,O,I,J) * ECOL(M,O,I,J)  *  &
                        ABS(vterm(O,J)-vterms) * pi * (RB+rs)**2.

              ! Collection kernel for interstitial scavenging
              ! if critcial radius is below 1 um
              ! upper radius limit rlim (= 200 nm)
                if ((rs.le.rlim).and.(rc.lt.1.e-6_dp)) then
                  if (rs.lt.2e-8_dp) then
                    KCOL(M,O,I,J) = KCOA(M,O,I,J)
                  else   ! 20 nm < rs < 400 nm
                    !KCOL(M,O,I,J) = KCOA(M,O,I,J)*(7.25e8*rs-12.5_dp)  ! factor 60
                    !KCOL(M,O,I,J) = KCOA(M,O,I,J)*(9.75e8*rs-17.5_dp)  ! factor 80
                    KCOL(M,O,I,J) = KCOA(M,O,I,J)*(1.23e9*rs-22.5_dp)  ! factor 100
                    !!!KCOL(M,O,I,J) = KCOA(M,O,I,J)*(1.48e9*rs-27.5_dp)  ! factor 120
                  endif
                endif

              ! No collection until critcial radius is below 1 um
                if (rc.ge.1.e-6_dp) then
                  KCOL(M,O,I,J) = 0._dp
                endif

              ! RB < rs
              ELSE
                ECOL(M,O,I,J) = 0._dp
                ECOA(M,O,I,J) = 0._dp
                KCOL(M,O,I,J) = 0._dp
              ENDIF

           !   if ((M==AI).and.(I==10).and.(rc.lt.1.e-6_dp)) then
           !   write(6,*) 'Kcol',O,J,rs/RB,ECOL(M,O,I,J),KCOA(M,O,I,J),KCOL(M,O,I,J)
           !   endif

            enddo   ! J
          enddo     ! O

        enddo       ! I
      enddo         ! M          


  end subroutine collection_kernel



    function gaussquad(n) result(r)
    !----------------------------------------------------------------------
    !
    !****
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Numerical integration / Gauss-Legendre Quadrature
    !      Works with gfortran but needs the option 
    !        -assume realloc_lhs
    !      when compiled with Intel Fortran.
    !
    !      interface
    !      ---------
    !
    !        input:
    !           n      number of Gaussian evaluation points
    !
    !      method
    !      ------
    !      The input values should be an function f to integrate, the bounds 
    !      of the integration interval a and b, and the number of gaussian 
    !      evaluation points n
    !
    !      external
    !      --------
    !      none
    !
    !      reference
    !      ---------
    !      https://rosettacode.org/wiki/Numerical_integration/
    !      Gauss-Legendre_Quadrature
    !
    !
    !------------------------------------------------------------------
    implicit none

    integer, intent(in)                       :: n

    integer                                   :: k, i,  iter
    real( dp)                                 :: x, f, df, dx
    real( dp), parameter                      :: pi = 4*atan(1._dp)

    real( dp), allocatable                    :: p0(:), p1(:), tmp(:)
    real( dp)                                 :: r(2, n)

    p0 = [1._dp]
    p1 = [1._dp, 0._dp]
 
    do k = 2, n
      tmp = ((2*k-1)*[p1,0._dp]-(k-1)*[0._dp, 0._dp,p0])/k
      p0 = p1
      p1 = tmp
    end do

    do i = 1, n
      x = cos(pi*(i-0.25_dp)/(n+0.5_dp))
      do iter = 1, 10
        f = p1(1)
        df = 0._dp
        do k = 2, size(p1)
          df = f + x*df
          f  = p1(k) + x * f
        end do
        dx =  f / df
        x = x - dx
        if (abs(dx)<10*epsilon(dx)) exit
      end do
      r(1,i) = x
      r(2,i) = 2/((1-x**2)*df**2)
    end do

  end function gaussquad

    
end module gde_coagulation
