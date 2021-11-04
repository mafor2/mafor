!
! <prog_fitaero.for - A component of FITAERO
!                 Fit utility for number size distributions>
!**********************************************************************!
!*
!**********************************************************************! 
!* 
!*    Copyright (C) 2011-2020  Matthias Steffen Karl
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
!*    1.  Karl, M., Gross, A., Pirjola, L., Leck, C., A new flexible
!*        multicomponent model for the study of aerosol dynamics
!*        in the marine boundary layer, Tellus B, 63(5),1001-1025,
!*        doi:10.1111/j.1600-0889.2011.00562.x, 2011.
!*    2.  Karl, M., Kukkonen, J., Keuken, M.P., Lutzenkirchen, S.,
!*        Pirjola, L., Hussein, T., Modelling and measurements of urban
!*        aerosol processes on the neighborhood scale in Rotterdam,
!*        Oslo and Helsinki, Atmos. Chem. Phys., 16,
!*        4817-4835, doi:10.5194/acp-16-4817-2016, 2016.
!*
!***********************************************************************! 
!*
!***********************************************************************
!***
!***      FITAERO
!***      Aerosol size distribution fit tool
!***      for MAFOR's inaero.dat
!***********************************************************************

      program prog_fitaero

!***********************************************************************
!***  This program for minimization of a function is based on SIMPLEX:
!***  
!***     https://oregonstate.edu/instruct/ch590/lessons/lesson10.htm
!***
!***  SIMPLEX is an optimization technique in which one finds the best 
!***  response in a multidimensional space, where the variables may be 
!***  interrelated.  Simple is computationally slow but is quick and 
!***  easy to implement.  The essence of the method is to observe the 
!***  response of a system, change variables intelligently, note the new
!***  response, change variables, etc.  One is searching for a minimum
!***  or a maximum in a response surface.  The response surface is in 
!***  n-dimensional space where the axes correspond to the variables of 
!***  the problem.  A simplex is a geometric figure with N+1 vertices 
!***  in an n-dimensional space.
!***
!***********************************************************************

!     Declarations of variables by using the MODULES feature:

      use module_fitaero_io
      use module_fitaero_exe

      implicit none

!***********************************************************************

!     Local Declarations:

! Input
      integer                     :: res1, res2
      integer                     :: resin1,resin2
      integer                     :: in
      character (len=256)         :: header1,header2
      logical                     :: leof

! Simplex
      integer                     :: nsim
      integer                     :: i,j,k,m,o
      integer                     :: ilo, ihi, nhi
      integer                     :: mstart
      real                        :: erfr, ex

! Size distribution
      integer, dimension(nm)      :: binlo
      integer, dimension(nm)      :: binup
      real, dimension(nm)         :: ntot
      real, dimension(nm)         :: vpt
      real, dimension(nm)         :: wamass
      real, dimension(nm)         :: shrmode
      integer                     :: nobsm
      real                        :: yhat
      real                        :: massbin
      real                        :: masswbin
      real                        :: massmode
      real                        :: masswmode
      real                        :: shrf
      real                        :: swfmode
      real                        :: aitmax
      real                        :: accmax


!***********************************************************************

!     Flow control

!     Get user-supplied meta information
      call get_user_input

!     Get file units
      funit_in_sizedis      = get_file_unit()

!     Allocate arrays
      if (.not. allocated(in_array_dp))   allocate(in_array_dp(n_obins))
      if (.not. allocated(in_array_dn))   allocate(in_array_dn(n_obins))
      if (.not. allocated(in_array_nuf))  allocate(in_array_nuf(nc))
      if (.not. allocated(in_array_aif))  allocate(in_array_aif(nc))
      if (.not. allocated(in_array_acf))  allocate(in_array_acf(nc))
      if (.not. allocated(in_array_cof))  allocate(in_array_cof(nc))
      if (.not. allocated(dataa))         allocate(dataa(n_obins,nv))
      if (.not. allocated(nest))          allocate(nest(nm,n_obins))
      if (.not. allocated(dense))         allocate(dense(nm,nc-1))


!     Initialize arrays with zeroes
      in_array_dp(:)    = 0.0
      in_array_dn(:)    = 0.0
      in_array_nuf(:)   = 0.0
      in_array_aif(:)   = 0.0
      in_array_acf(:)   = 0.0
      in_array_cof(:)   = 0.0
      dataa(:,:)        = 0.0
      dense(:,:)        = 0.0
      wamass(:)         = 0.0

      kount=0

      print *,'FITDIS after get_user_input'

!     Load Initial Distribution Data
!     bin-number, Dp [nm], dN [#/cm^3]

      open(funit_in_sizedis, file=fname_in_sizedis,access='sequential',&
           form="formatted",iostat=res1)
      read(funit_in_sizedis, fmt='(A)', iostat=resin1) header1
      call read_csv_file(funit_in_sizedis, n_obins, in_array_dp)
      call read_csv_file(funit_in_sizedis, n_obins, in_array_dn)


!     Load Mass Fraction data
!     H2SO4 ORG NH4 NO3 MSAp SALT BPAB BC DUST Total
      open(funit_in_massfra,file=fname_in_massfra,access='sequential',&
           form="formatted",iostat=res2)
      read(funit_in_massfra, fmt='(A)', iostat=resin2) header2
      call read_csv_file(funit_in_massfra, nc, in_array_nuf)
      call read_csv_file(funit_in_massfra, nc, in_array_aif)
      call read_csv_file(funit_in_massfra, nc, in_array_acf)
      call read_csv_file(funit_in_massfra, nc, in_array_cof)


      if (fe_log) then
        write(funit_log,*)
        write(funit_log,'(1X,A58)')  'FITDIS after input' 
        write(funit_log,*)
      endif

!     Check values in data matrix
      do i = 1, n_obins
         if (in_array_dn(i)==0) then
           call stopit_fitaero("Observed dN value must not be zero")
         endif
      end do

!     Fill data matrix (x,y)
      do i = 1, n_obins
         dataa(i,1) = in_array_dp(i)
         dataa(i,2) = in_array_dn(i)
      end do

!     Fill mass frac matrix
      do o = 1, nc-1
         dense(1,o) = in_array_nuf(o)
         dense(2,o) = in_array_aif(o)
         dense(3,o) = in_array_acf(o)
         dense(4,o) = in_array_cof(o)
      end do

!     Fill parameter estimates
!     Parameters' upper limit
!     Parameters' lower limit
!     Use FINE set for Dpmax < 2.0E-6 m
!     Use COARSE+FINE set for 2.0E-6 >= Dpmax < 1.E-5 m
!     Use ARCTIC set for Dpmax = 1.E-5 m 
      do m = 1, nm
         if (dpmax<2.0E-6) then
           beta(m,1)    = es2_gmd(m)
           beta(m,2)    = es2_sigma(m)
           beta(m,3)    = es2_mm(m)
           betamax(m,1) = max2_gmd(m)
           betamax(m,2) = max2_sigma(m)
           betamax(m,3) = max2_mm(m)
           betamin(m,1) = min2_gmd(m)
           betamin(m,2) = min2_sigma(m)
           betamin(m,3) = min2_mm(m)
         else if ( (dpmax>=2.0E-6).and.(dpmax<1.0E-5) ) then
           beta(m,1)    = es1_gmd(m)
           beta(m,2)    = es1_sigma(m)
           beta(m,3)    = es1_mm(m)
           betamax(m,1) = max1_gmd(m)
           betamax(m,2) = max1_sigma(m)
           betamax(m,3) = max1_mm(m)
           betamin(m,1) = min1_gmd(m)
           betamin(m,2) = min1_sigma(m)
           betamin(m,3) = min1_mm(m)
         else
           beta(m,1)    = es3_gmd(m)
           beta(m,2)    = es3_sigma(m)
           beta(m,3)    = es3_mm(m)
           betamax(m,1) = max3_gmd(m)
           betamax(m,2) = max3_sigma(m)
           betamax(m,3) = max3_mm(m)
           betamin(m,1) = min3_gmd(m)
           betamin(m,2) = min3_sigma(m)
           betamin(m,3) = min3_mm(m)
         endif
      end do

      print *,'betaes ',beta(:,1)
      print *,'betamax',betamax(:,1)      
      print *,'betamin',betamin(:,1)


!     Nucleation mode yes/no
      if ( dataa(1,1).lt.9.0 ) then
         mstart=1
      else
         mstart=2
         binlo(2)=1
      endif

!     Find lower and upper bin of each mode
      call findmodes(dataa,n_obins,mstart,binlo,binup)

      print *,'bin low',binlo(:)
      print *,'bin  up',binup(:)


!     Update the lower and upper limit for gmd
      do m = 2, 3
        betamax(m,1) = dataa(binup(m),1)
        betamin(m,1) = dataa(binlo(m),1)
      enddo

!     Adjust Aitken mode limits
      if ( (dpmax>=2.0E-6).and.(dpmax<1.0E-5) ) then
         betamin(2,1) = min1_gmd(2)
      endif
      betamax(2,1) =betamax(2,1)*0.91


!     Set gmd limits for COARSE mode
      if (dpmax<2.0E-6) then
        betamin(4,1) = min2_gmd(4)
      else if ( (dpmax>=2.0E-6).and.(dpmax<1.0E-5) ) then
        betamin(4,1) = min1_gmd(4)         
      else
        betamin(4,1) = min3_gmd(4)
      endif
      betamax(4,1) = dataa(binup(4),1)


!     Find Aitken and Acc mode maximum
      aitmax = maxval( dataa(binlo(2):binup(2),2) )
      accmax = maxval( dataa(binlo(3):binup(3),2) )
      
      do i = 1, n_obins
        if (dataa(i,2).eq.aitmax) then
           beta(2,1) = dataa(i,1)
        endif
        if (dataa(i,2).eq.accmax) then
           beta(3,1) = dataa(i,1)
        endif
      enddo


!     Estimate of mass based on observed dN and GMD_es
!     MASS[ng/m3] = N[1/cm3]*VPT(gmd)[m3]*CONV[ng/kg]*DENS[kg/m3]*1.E6
      do m = 1, nm
        !print *,'gmd',m,beta(m,1)
        ntot(m)    = sum(dataa(binlo(m):binup(m),2))
        vpt(m)     = (pi/6.)*(beta(m,1)*1.e-9)**3.0
        beta(m,3)  = ntot(m)*vpt(m)*conv*dens(m)*1.e6
        shrmode(m) = 1.0
        print *,'es',m,ntot(m),vpt(m),beta(m,3)
      enddo

      betamin(:,3)=betamin(:,3)*0.08

!     Update upper mass limit
      if (dpmax<1.0E-5)  then
        do m = mstart, nm
          betamax(m,3) = beta(m,3)*300.0
        enddo
      endif

!     Update lower mass limit (AI Background)
      if ( (dpmax>=2.0E-6).and.(dpmax<1.0E-5) ) then
        if(beta(2,3)<50.0) then
          if (dense(2,SU)>0.4) then 
            betamin(2,3) =  5.0  ! ug/m^3 (marine)
          else
            betamin(2,3) = 20.0  ! ug/m^3 (coastal)
          endif
        else if(beta(2,3)<100.0) then   !@ 82 ug/m3
          betamin(2,3) = 280.0  ! ug/m^3
          betamin(2,3) =  70.0  ! ug/m^3 (marine)
        else if(beta(2,3)<200.0) then
          if (dense(2,SU)>0.4) then 
            betamin(2,3) = 70.0  ! ug/m^3 (marine)
          else
            betamin(2,3) = 400.0  ! ug/m^3 (coastal)
          endif
        else if(beta(2,3)<300.0) then
          betamin(2,3) = 200.0  ! ug/m^3 (coastal)
        endif
      endif


!     Generic Nucleation Mode (mstart=2)
!     Add small mass in nucleation mode if zero
      if (mstart==2) then
         beta(1,1) =  9.5  ! nm
         beta(1,2) =  1.30
         beta(1,3) =  0.45 * dataa(1,2)/1200.0
      endif

!     Adjustment of Exhaust mode widths
      if (beta(2,3)>10000.0) then
         betamax(1,1) = 7.0 ! nm
         betamax(1,2) = 1.35 
         betamin(2,2) = 1.62
         betamax(2,2) = 1.85
         betamin(3,2) = 1.45
         betamax(3,2) = 1.65
      endif


! ** Check here the parameter overview
      print *,'gmd min',betamin(:,1)
      print *,'gmd max',betamax(:,1)
      print *,'gmd est',beta(:,1)

      print *,'sig min',betamin(:,2)
      print *,'sig max',betamax(:,2)
      print *,'sig est',beta(:,2)

      print *,'mm min',betamin(:,3)
      print *,'mm est',beta(:,3)
      print *,'mm max',betamax(:,3)

    !! stop

! **************************
! *** BEGIN THE SIMPLEX  ***
! **************************

      do m = mstart, nm

       print *,'new mode',m
       print *,'bins in mode',binup(m)-binlo(m)+1

       nobsm = binup(m)-binlo(m)+1
       nsim = 1
       erf(:)=0.0


       print *,'lo up1',m,binlo(m),binup(m)
       print *,'errf1',m,(beta(m,j),j=1,np)
       call errfunc(beta,dataa,dense,nv,binlo,binup,m,kount,erf(1))

          
       if (fe_log) then
         write(funit_log,*)
         write(funit_log,'(1X,A58,G12.4,A6,I2)') 'start error valu ', &
                                           erf(1),' mode ',m
         write(funit_log,*)
       endif



! *** Initialize the simplex
       kount = 0

       do 22 j=1,np
          p(m,1,j) = beta(m,j)
 22    end do

       do 28 i=2,np1

         do 26 j=1,np
            p(m,i,j) = beta(m,j)
 26      end do

         p(m,i,i-1) = 1.1 * beta(m,i-1)
         if ( abs(beta(m,i-1)).lt.1.0e-12 ) then
            p(m,i,i-1) = 0.001
         endif

 28    continue


       if (fe_log) then
         write(funit_log,*)
         write(funit_log,'(1X,A58,i2)')  'FITDIS after init, mode ',m 
         write(funit_log,*)
       endif




! *** Start walking the simplex

! find p_low and p_high  / best=p_low / worst=p_high
 31    ilo=1
       ihi=1

       do 34 i=1,np1

         do 32 j=1,np

           beta(m,j)=p(m,i,j)

 32      end do

! *** Constrain betas with upper & lower limits
         do j=1,np
            beta(m,j)=min(beta(m,j),betamax(m,j))
            beta(m,j)=max(beta(m,j),betamin(m,j))
         enddo


         !print *,'lo up2',m,binlo(m),binup(m)
         !print *,'errf2',m,(beta(m,j),j=1,np)
         call errfunc(beta,dataa,dense,nv,binlo,binup,m,kount,erf(i))

         if (erf(i).lt.erf(ilo)) ilo=i
         if (erf(i).gt.erf(ihi)) ihi=i

 34    continue

       if (fe_log) then
         write(funit_log,*)
         write(funit_log,'(1X,A58)')  'FITDIS initial simplex'
         write(funit_log,*)
         write(funit_log,'(1X,A)') '********************************* '

         do 40 k=1,np1
           write(funit_log,39)k,erf(k), (p(m,k,j),j=1,np)
 40      continue

         do k=1,np1
          print *,'ini simplex',(p(m,k,j),j=1,np)
         enddo
       endif


! find p_nhi the next highest next=p_nhi     
 41    nhi=ilo

       do 43 i=1,np1

         if (erf(i).ge.erf(nhi).and.i.ne.ihi) nhi=i

 43    continue

! compute the centroid
       do 46 j=1,np
         c(m,j)=(-1)*p(m,ihi,j)
         do 44 i=1,np1
           c(m,j)=c(m,j)+p(m,i,j)
 44      continue
         c(m,j)=c(m,j)/np
 46    continue
 
 51    continue

! print current best vertex

       if (fe_log) then
         write(funit_log,'(1X,A)') '********************************* '
         write(funit_log,*)
         write(funit_log,'(1X,A58)')  'FITDIS current best vertex'
         write(funit_log,*)
         write(funit_log,53)kount,nsim
         write(funit_log,54)(p(m,ilo,j),j=1,np)
         !write(funit_log,'(1X,A58,G12.4)') 'error function', erf(ilo)
       endif
 
 
! stopping criterion

       if (kount.gt.maxcnt) then
         print *,'!!!max. number simulations reached for mode ',m
         !call stopit_fitaero('Max. number of simulations reached')
         if (fe_log) then
         write(funit_log,'(1X,A)') '********************************* '
         write(funit_log,*)
         write(funit_log,'(1X,A58)')  'FITDIS max.number of simulations'
         write(funit_log,*)
         endif
       endif

! for the next condition, erf(ilo) must not be zero
       erf(ilo) = max(erf(ilo),1.e-18)

       if(( (abs(erf(ilo)-erf(ihi))/erf(ilo)).lt.errmin) &
            .or.(kount.gt.maxcnt)) then
      
        ! End of simplex

 56     print *,'  ===> error criterion satisfied'

        if (fe_log) then
          write(funit_log,54)(p(m,ilo,j),j=1,np)
          write(funit_log,*)
          write(funit_log,'(1X,A58)') 'FITDIS error criterion fulfilled'
          write(funit_log,*)
        endif

        print *,'Normal termination of simplex in mode ',m


        ! Write output
        do j=1,np
           beta(m,j)=p(m,ilo,j)
        enddo
        shrf=1.0
        massmode=0.0
        masswmode=0.0
        do i= binlo(m),binup(m)

           call initaero(beta,dataa,dense,shrf,binlo,binup,m,i,  &
                        yhat,massbin,masswbin)

           massmode=massmode+massbin
           masswmode=masswmode+masswbin

           nest(m,i)=yhat
           !print *,'dNest ',m,i,dataa(i,1),dataa(i,2),nest(m,i)
        enddo
        ! fraction of water in mode
        wamass(m)=(masswmode-massmode)/masswmode


        ! output the shrinking factor
        swfmode=(masswmode/(max(massmode,1.e-20)))**(0.3333)
        swfmode=min(swfmode,3.0)
        swfmode=max(swfmode,1.0)
        shrmode(m)=1./swfmode

        print *,'shr',shrmode(m)

        do j=1,np
           print *,'beta before output m', m,beta(m,j)
        enddo

        call output_fitaero(beta,dataa,dense,ntot,wamass,shrmode,  &
                            binlo,binup,m,mstart,nest)

       if (m.le.nm) then
          go to 140
        else
          call stopit_fitaero('End of simplex')
        endif

       endif

! *** Reflection of the simplex   

 61    do 62 j=1,np
         r(m,j)=1.9985*c(m,j)-0.9985*p(m,ihi,j)
 62    continue

! *** Constrain betas with upper & lower limits
       do j=1,np
          r(m,j)=min(r(m,j),betamax(m,j))
          r(m,j)=max(r(m,j),betamin(m,j))
       enddo

       !print *,'lo up3',m,binlo(m),binup(m)
       !print *,'errf3',m,(r(m,j),j=1,np) 
       call errfunc(r,dataa,dense,nv,binlo,binup,m,kount,erfr)


       if (fe_log) then
         write(funit_log,65)erfr,(r(m,j),j=1,np)
       endif

! reflect again if successful

       if(erfr.lt.erf(ilo)) go to 91    ! expand simplex

       if(erfr.ge.erf(ihi)) go to 122   ! contract simplex

! otherwise replace worst vertex with new one

 79    do 80 j=1,np     
         p(m,ihi,j)=r(m,j)
 80    continue

       ! increase number of simulations
       nsim=nsim+1

       erf(ihi)=erfr

       if(erfr.gt.erf(nhi)) go to 51

       ihi=nhi

       go to 41


! *** Expand the simplex  

 91    ilo=ihi

       ihi=nhi

       do 93 j=1,np
          beta(m,j) = 1.95*r(m,j) - 0.95*c(m,j)
 93    continue


! *** Constrain betas with upper & lower limits
       do j=1,np
          beta(m,j)=min(beta(m,j),betamax(m,j))
          beta(m,j)=max(beta(m,j),betamin(m,j))
       enddo

       !print *,'lo up4',m,binlo(m),binup(m)
       !print *,'errf4',m,(beta(m,j),j=1,np) 
       call errfunc(beta,dataa,dense,nv,binlo,binup,m,kount,ex)


       if(ex.lt.erfr) go to 104


! r is better than beta
       do 99 j=1,np
          p(m,ilo,j) = r(m,j)
 99    continue

      ! increase number of simulations
       nsim = nsim + 1

       erf(ilo)=erfr

       go to 110



! beta is better than r
 104   do 105 j=1,np     
          p(m,ilo,j) = beta(m,j)
 105   continue

       if (fe_log) then
         write(funit_log,106)ex,(beta(m,j),j=1,np)
         write(funit_log,*)
         write(funit_log,'(1X,A58)')  'FITDIS expansion vortex done.'
         write(funit_log,*)
       endif


       ! increase number of simulations
       nsim=nsim+1

       erf(ilo)=ex

 110   continue

       go to 41


! *** Contract the simplex  

 122   do 123 j=1,np     
         r(m,j)=0.5015*c(m,j)+0.4985*p(m,ihi,j)
 123   continue

! *** Constrain betas with upper & lower limits
       do j=1,np
          r(m,j)=min(r(m,j),betamax(m,j))
          r(m,j)=max(r(m,j),betamin(m,j))
       enddo

       !print *,'lo up5',m,binlo(m),binup(m)
       !print *,'errf5',m,(r(m,j),j=1,np) 
       call errfunc(r,dataa,dense,nv,binlo,binup,m,kount,erfr)

       if (fe_log) then
         write(funit_log,124) erfr, (r(m,j),j=1,np)
         write(funit_log,*)
         write(funit_log,'(1X,A58)')  'FITDIS contraction vortex done.'
         write(funit_log,*)
       endif

       if(erfr.lt.erf(ilo)) go to 91

       if(erfr.lt.erf(ihi)) go to 79

! *** Scaling
       if (scalf.ne.0.0) then

 137     do 138 i=1,np1
           do 139 j=1,np
              p(m,i,j)=p(m,i,j)+scalf*(p(m,ilo,j)-p(m,i,j))
 139       continue
 138     continue

         print*,'hello start walking again'
         go to 31   ! start walking again

       endif

! *** End of simplex algorithm


 140   print *,'end mode',m

!     End of loop over aerosol modes
      end do


! NOW WE SHOULD WRITE THE OUTPUT FILE

!     Close all i/o files
      close(funit_in_sizedis)
      close(funit_in_massfra)
      close(funit_out_sizedis)

!     Deallocate variables 
      deallocate(in_array_dp)
      deallocate(in_array_dn)
      deallocate(in_array_nuf)
      deallocate(in_array_aif)
      deallocate(in_array_acf)
      deallocate(in_array_cof)
      deallocate(dataa)
      deallocate(nest)
      deallocate(dense)

      if (fe_log) then
        write(funit_log,*)
        write(funit_log,'(1X,A58)')  'FITDIS finished'
        write(funit_log,*)
      endif

!     Close the Log-file
      close(funit_log) 

!     Format statements
 39   format(3x,' vertex',i2,' eror and paramerters:',5f12.3)
 53   format(' after ', i4, ' error evaluations and ',i4,' simplexes')
 54   format('  parameterm estimates:    ',4f12.3)
 65   format(' reflection vertex', 4f12.3)
 106  format(' expansion  vertex', 4f12.3)
 124  format(' contracton vertex', 4f12.3)

      end program prog_fitaero
