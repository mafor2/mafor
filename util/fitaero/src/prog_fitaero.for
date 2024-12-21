!
! <prog_fitaero.for - A component of FITAERO
!                 Fit utility for number size distributions>
!**********************************************************************!
!*
!**********************************************************************! 
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
!***      Version: 1.3
!***
!***      FITAERO treats three sets of aerosols
!***      dpmax <  2.0E-6 [m]            (FINE)
!***      2.0E-6 <= dpmax < 1.0E-5  [m]  (FINE + COARSE)
!***      dpmax = 1.0E-5  [m]            (PM10)
!***
!***      FINE subgroups:
!***          marine  background [beta(3,3) <   500.0]
!***          coastal background [beta(3,3) <  2000.0]
!***          urban   background [beta(3,3) < 10000.0]
!***          bioburn background [beta(3,3) < 25000.0]
!***          exhaust/combustion [beta(3,3) > 25000.0]
!***
!***      FINE+COARSE subgroups:
!***          marine  background [beta(2,3) <   100.0]
!***          coastal background [beta(2,3) <   200.0]
!***          urban   background [beta(2,3) < 10000.0]
!***          exhaust/combustion [beta(2,3) > 10000.0]
!***
!***      PM10 subgroups:
!***          marine  background [beta(2,3) <   200.0]
!***          urban   background [beta(2,3) < 10000.0]
!***          exhaust/combustion [beta(2,3) > 10000.0]
!***
!***
!***      How to add a new Subgroup
!***      prog_fitaero:       L375-600  gmd, sigma, mass estimate
!***      module_fitaero_exe: L333-451  shrmode --> mass change 
!***      output_fitaero:     L173-266  gmd, sigma after fitting
!***
!***      How to calculate dlogDp of observed Number Size Distr
!***          dlogDp = LOG10 ( Dp(i+1) / Dp(i) )
!***
!***********************************************************************
!*
!***********************************************************************
!***
!***      REVISION HISTORY
!***
!***  15 Jun 2023 M. Karl L763 close(funit_out_inaero)
!***  17 Jun 2023 M. Karl      maxcnt = 20 error evaluations (m=2)
!***  17 Jun 2023 M. Karl      improved mass estimate calculation
!***  17 Jun 2023 M. Karl      output_fitaero subtract water mass (m>2)
!***  29 Dec 2023 M. Karl      output_fitaero subtract water only m=4
!***  31 Dec 2023 M. Karl      module_fitaero_exe same sigma estimates
!***                           for all sets, as in FITAERO_params.xls
!***  03 Jan 2024 M. Karl      FINE+COARSE subgroups coastal & marine
!***  04 Jan 2024 M. Karl      PM10 subgroups marine & urban
!***  13 Jan 2024 M. Karl      FINE subgroup exhaust
!***  20 Jan 2024 M. Karl      reset mode diameter lower limits after
!***                           call to findmodes
!***  06 Dec 2024 M. Karl      extend FITAERO to 5 mode-distribution
!***  19 Dec 2024 M. Karl      recalculate output sizedis with lognormal
!***                           mass concentration distribution
!***
!*** TO DO: remove shrmode adjustments in module_fitaero_exe.for
!***        unless necessary to distinguish cases
!***
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
      integer                     :: i,j,k,m,o,b,n
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
      real                        :: mass
      real                        :: numb
      real                        :: yhat
      real                        :: massbin
      real                        :: masswbin
      real                        :: massmode
      real                        :: masswmode
      real                        :: shrf
      real                        :: swfmode
      real                        :: napmax
      real                        :: aitmax
      real                        :: accmax

! Output size distribution (inaero.dat)
      real, dimension(nm)         :: gmdout
      real, dimension(nm)         :: sigout
      real                        :: massout
      real                        :: volcout


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
      if (.not. allocated(in_array_naf))  allocate(in_array_naf(nc))
      if (.not. allocated(in_array_aif))  allocate(in_array_aif(nc))
      if (.not. allocated(in_array_acf))  allocate(in_array_acf(nc))
      if (.not. allocated(in_array_cof))  allocate(in_array_cof(nc))
      if (.not. allocated(dataa))         allocate(dataa(n_obins,nv))
      if (.not. allocated(nest))          allocate(nest(nm,n_obins))
      if (.not. allocated(numout))        allocate(numout(nm,n_obins))
      if (.not. allocated(dense))         allocate(dense(nm,nc-1))


!     Initialize arrays with zeroes
      in_array_dp(:)    = 0.0
      in_array_dn(:)    = 0.0
      in_array_nuf(:)   = 0.0
      in_array_naf(:)   = 0.0
      in_array_aif(:)   = 0.0
      in_array_acf(:)   = 0.0
      in_array_cof(:)   = 0.0
      dataa(:,:)        = 0.0
      dense(:,:)        = 0.0
      wamass(:)         = 0.0
      gmdout(:)         = 0.0
      sigout(:)         = 0.0

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
      call read_csv_file(funit_in_massfra, nc, in_array_naf)
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
         dense(2,o) = in_array_naf(o)
         dense(3,o) = in_array_aif(o)
         dense(4,o) = in_array_acf(o)
         dense(5,o) = in_array_cof(o)
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

      print *,'size ranges (nm)' 
      print *,'betamin',betamin(:,1)
      print *,'betamax',betamax(:,1)
      print *,'betaes ',beta(:,1)


!     Nucleation mode yes/no
!     Include fitted NU mode if lowest
!     measured particle diameter is <7 nm
      if ( dataa(1,1).lt.7.0 ) then
         mstart=1
      else
         mstart=2
         binlo(2)=1
      endif

!     Find lower and upper bin of each mode
      call findmodes(dataa,n_obins,mstart,binlo,binup)

      print *,'start mode',mstart
      print *,'bin low',binlo(:)
      print *,'bin  up',binup(:)
      !stop

!     Update the lower and upper limit for gmd
      do m = NA, AS
        betamax(m,1) = dataa(binup(m),1)
        betamin(m,1) = dataa(binlo(m),1)
      enddo
      betamax(CS,1) = dataa(binup(CS),1)


!     Reset Mode diameter lower limits
      !NAP,AIT
      if ( (dpmax>=2.0E-6).and.(dpmax<1.0E-5) ) then
         betamin(NA,1) = min1_gmd(2)
         betamin(AI,1) = min1_gmd(3)
      endif
      !ACC,COA
      if (dpmax<2.0E-6) then
        betamin(AS,1) = min2_gmd(4)
        betamin(CS,1) = min2_gmd(5)
      else if ( (dpmax>=2.0E-6).and.(dpmax<1.0E-5) ) then
        betamin(AS,1) = min1_gmd(4)  
        betamin(CS,1) = min1_gmd(5)
      else
        betamin(AS,1) = min3_gmd(4)
        betamin(CS,1) = min3_gmd(5)
      endif


!     Find Nanoparticle, Aitken and Acc mode maximum
      napmax = maxval( dataa(binlo(NA):binup(NA),2) )
      aitmax = maxval( dataa(binlo(AI):binup(AI),2) )
      accmax = maxval( dataa(binlo(AS):binup(AS),2) )

     
      do i = 1, n_obins
        if (dataa(i,2).eq.napmax) then
           beta(NA,1) = dataa(i,1)
        endif
        if (dataa(i,2).eq.aitmax) then
           beta(AI,1) = dataa(i,1)
        endif
        if (dataa(i,2).eq.accmax) then
           beta(AS,1) = dataa(i,1)
        endif
      enddo

      print *,'beta GMD',beta(:,1)

!     Estimate of mass based on observed dN and GMD_es
!     MASS[ng/m3] = N[1/cm3]*VPT(gmd)[m3]*CONV[ng/kg]*DENS[kg/m3]*1.E6
      do m = 1, nm
        ntot(m)    = 0.0
        beta(m,3)  = 0.0
        mass       = 0.0
        numb       = 0.0
        do b = binlo(m), binup(m)
          numb       = (dataa(b,2))
          vpt(m)     = (pi/6.)*((dataa(b,1))*1.e-9)**3.0
          mass       = numb*vpt(m)*conv*dens(m)*1.e6
          beta(m,3)  = beta(m,3) + mass
          ntot(m)    = ntot(m)   + numb
        enddo
        print *,'mass es',m,beta(m,3)
        shrmode(m) = 1.0
      enddo


!     Update upper mass limits and mass estimates

! FINE set
      if (dpmax<2.0E-6) then

         print *,'FINE set'

         print *,'mstart ',mstart

         print *,'Mtot(AS)',beta(AS,3)

! ** Division in exhaust / urban backgr / coastal backgr

        if (beta(AS,3)<500.0) then
        !bkgr2_test

           print *,'MARINE'
           !NAP
           beta(NA,1)     = beta(NA,1)   *0.85
           betamin(NA,1)  = betamin(NA,1)*0.85
           betamax(NA,1)  = betamax(NA,1)*0.85
           !AIT
           beta(AI,1)     = beta(AI,1)   *1.20
           betamin(AI,1)  = betamin(AI,1)*1.20           
           betamax(AI,1)  = betamax(AI,1)*1.20
           beta(AI,3)     = beta(AI,3)   *0.73
           !ACC
           beta(AS,1)     = beta(AS,1)   *1.2
           betamin(AS,1)  = betamin(AS,1)*1.2
           betamax(AS,1)  = betamax(AS,1)*1.2

           ! change mode width
           betamin(AI,2)  = 1.87
           betamax(AI,2)  = 1.99
           !update mass min & max
           betamin(NU:NA,3)=betamin(NU:NA,3)/25.
           betamin(AI,3)  = betamin(AI,3)   *0.15
           betamin(AS:CS,3)=betamin(AS:CS,3)*0.3
           betamax(:,3)   = betamax(:,3)/100.

        else if (beta(AS,3)<3000.0) then
        !bktr1_test

           print *,'COASTAL'
           !NAP
           beta(NA,1)     = beta(NA,1)   *0.90
           betamin(NA,1)  = betamin(NA,1)*0.90
           betamax(NA,1)  = betamax(NA,1)*0.90
           beta(NA,3)     = beta(NA,3)   *0.50
           !AIT
           beta(AI,1)     = beta(AI,1)   *0.90
           betamin(AI,1)  = betamin(AI,1)*0.90
           betamax(AI,1)  = betamax(AI,1)*0.90
           !ACC
           beta(AS,1)     = beta(AS,1)   *1.40
           betamin(AS,1)  = betamin(AS,1)*1.40
           betamax(AS,1)  = betamax(AS,1)*1.40

           ! change mode width
           betamin(NA,2)  = 1.55
           betamax(NA,2)  = 1.65
           betamin(AI,2)  = 1.82
           betamax(AI,2)  = 1.98
           betamin(AS,2)  = 1.85
           betamax(AS,2)  = 1.95
           !update mass min & max
           betamin(NA,3)  = betamin(NA,3)*0.50
           betamax(NA,3)  = betamax(NA,3)*0.50
           beta(AI:CS,3)  = beta(AI:CS,3)*10.0
           betamin(AI,3)  = beta(AI,3)   *0.17
           betamin(AS,3)  = beta(AS,3)   *0.16
           betamin(CS,3)  = beta(CS,3)   *0.25
           betamax(CS,1)  = max2_gmd(CS) *0.5

        else if ( beta(AS,3) < 11000.0) then
        !traff1_test: do not change parameters
        !init1_test:  do not change here but in output_fitaero

           print *,'URBAN'
           !NAP
           beta(NA,1)     = beta(NA,1)   *0.94
           betamin(NA,1)  = betamin(NA,1)*0.94
           betamax(NA,1)  = betamax(NA,1)*0.94
           beta(NA,3)     = beta(NA,3)   *2.5
           !AIT
           beta(AI,1)     = beta(AI,1)   *1.04
           betamin(AI,1)  = betamin(AI,1)*1.04
           betamax(AI,1)  = betamax(AI,1)*1.04
           !ACC
           beta(AS,1)     = beta(AS,1)   *1.3
           betamin(AS,1)  = betamin(AS,1)*1.3
           betamax(AS,1)  = betamax(AS,1)*1.3
           !COA
           betamin(CS,1)  = betamin(CS,1)*1.2

           ! change mode width
           betamin(NA,2)  = 1.60
           betamax(NA,2)  = 1.70
           betamin(AI,2)  = 1.82
           betamax(AI,2)  = 1.90
           !update mass min & max
           betamin(NA,3)  = betamin(NA,3)*2.0
           betamax(NA,3)  = betamax(NA,3)*2.0
           beta(AI:CS,3)  = beta(AI:CS,3)*3.0
           betamax(AI:CS,3) = beta(AI:CS,3)*4.0
           betamin(AI,3)  = beta(AI,3)   *0.77
           betamin(AS,3)  = beta(AS,3)   *0.49
           betamin(CS,3)  = beta(CS,3)   *0.45

        else if ( beta(AS,3) < 35000.0) then
        !aces_test

           print *,'BIOBURN'
           !NAP
           beta(NA,1)     = beta(NA,1)   *1.05
           betamin(NA,1)  = betamin(NA,1)*1.05
           betamax(NA,1)  = betamax(NA,1)*1.05
           !AIT
           beta(AI,1)     = beta(AI,1)   *0.80
           betamin(AI,1)  = betamin(AI,1)*0.80
           betamax(AI,1)  = betamax(AI,1)*0.80
           beta(AI,3)     = beta(AI,3)   *5.50
           !ACC
           beta(AS,1)     = beta(AS,1)   *1.06
           betamin(AS,1)  = betamin(AS,1)*1.06
           betamax(AS,1)  = betamax(AS,1)*1.06
           beta(AS,3)     = beta(AS,3)   *2.00

           ! change mode width
           betamin(NA,2)  = 1.45
           betamax(NA,2)  = 1.55
           betamin(AS,2)  = 1.60
           betamax(AS,2)  = 1.66
           !update mass min & max
           betamin(NA,3)  = betamin(NA,3)*0.75
           betamin(AI,3)  = betamin(AI,3)*50.
           betamin(AS:CS,3)  = betamin(AS:CS,3)*290.

        else
        !stena_test

           print *,'EXHAUST'
           !NUC
           beta(NU,1)     = beta(NU,1)   *1.13
           betamax(NU,1)  = betamax(NU,1)*1.13
           betamin(NU,1)  = betamin(NU,1)*1.13
           !NAP
           beta(NA,1)     = beta(NA,1)   *0.98
           betamax(NA,1)  = betamax(NA,1)*0.98
           betamin(NA,1)  = betamin(NA,1)*0.98
           !AIT
           beta(AI,3)     = beta(AI,3)   *2.1
           !ACC
           beta(AS,1)     = beta(AS,1)   *1.50
           betamin(AS,1)  = betamin(AS,1)*1.50
           betamax(AS,1)  = betamax(AS,1)*1.50
           !COA
           beta(CS,1)     = beta(CS,1)   *1.40
           betamin(CS,1)  = betamin(CS,1)*1.40
           betamax(CS,1)  = betamax(CS,1)*1.40
           beta(CS,3)     = beta(CS,3)   *0.21

           ! change mode width
           betamax(NU,2)  = 1.38
           betamin(NA,2)  = 1.43
           betamax(NA,2)  = 1.58
           betamin(AI,2)  = 1.62
           betamax(AI,2)  = 1.80
           betamin(AS,2)  = 1.58
           betamax(AS,2)  = 1.68
           betamin(CS,2)  = 1.70
           betamax(CS,2)  = 1.90
           !update mass min & max
           betamax(NA:CS,3)  = betamax(NA:CS,3)*300.
           betamin(NA,3)  = beta(NA,3)*0.62
           betamin(AI,3)  = beta(AI,3)*0.65
           betamin(AS,3)  = beta(AS,3)*1.08
           betamin(CS,3)  = beta(CS,3)*0.45

        endif


! 'FINE + COARSE set'
      else if ( (dpmax>=2.0E-6).and.(dpmax<1.0E-5) ) then

        print *,'FINE + COARSE set'

        print *,'mstart ',mstart
        print *,'Mtot(AI) ',beta(AI,3)

        ! update MASS min & max
        betamin(:,3) = betamin(:,3)*0.08

        do m = mstart, nm
           betamax(m,3) = beta(m,3)*300.0
        enddo

        if(beta(AI,3)<600.0) then
        !amar2_test
        ! ("coastal marine")
           print *,'MARINE'
           !NAP
           beta(NA,1)     = beta(NA,1)   *1.17
           betamax(NA,1)  = betamax(NA,1)*1.17
           betamin(NA,1)  = betamin(NA,1)*1.17
           beta(NA,3)     = beta(NA,3)*4.7
           !AIT
           beta(AI,1)     = beta(AI,1)   *1.07
           betamin(AI,1)  = betamin(AI,1)*1.07
           betamax(AI,1)  = betamax(AI,1)*1.07
           beta(AI,3)     = beta(AI,3)*2.1
           !ACC
           beta(AS,1)     = beta(AS,1)   *0.9
           betamin(AS,1)  = betamin(AS,1)*0.9
           betamax(AS,1)  = betamax(AS,1)*0.9
           beta(AS,3)     = beta(AS,3)*0.6
           !COA
           beta(CS,1)     = beta(CS,1)   *1.6
           betamax(CS,1)  = betamax(CS,1)*1.6
           betamin(CS,1)  = betamin(CS,1)*1.6
           beta(CS,3)     = beta(CS,3)*9.0

           ! change mode width
           betamin(AI,2)  = 1.75
           betamax(AI,2)  = 1.90
           !update mass min & max
           betamin(:,3)  = beta(:,3)*0.90

        else if (beta(AI,3)<2000.0) then
        !acont_test
        ! ("coastal continental")
           print *,'COASTAL'
           !NAP
           beta(NA,1)     = beta(NA,1)   *1.2
           betamax(NA,1)  = betamax(NA,1)*1.2
           betamin(NA,1)  = betamin(NA,1)*1.2
           beta(NA,3)     = beta(NA,3)*7.0
           !AIT
           beta(AI,1)     = beta(AI,1)   *2.0
           betamin(AI,1)  = betamin(AI,1)*2.0
           betamax(AI,1)  = betamax(AI,1)*2.0
           beta(AI,3)     = beta(AI,3)*9.4
           !ACC
           beta(AS,3)     = beta(AS,3)*1.0
           !COA
           beta(CS,1)     = beta(CS,1)   *1.6
           betamax(CS,1)  = betamax(CS,1)*1.6
           betamin(CS,1)  = betamin(CS,1)*1.6
           beta(CS,3)     = beta(CS,3)*9.0

           ! change mode width
           betamin(AI,2)  = 1.90
           betamax(AI,2)  = 2.15
           !update mass min & max
           betamin(:,3)  = beta(:,3)*0.90

        else if (beta(AI,3)<10000.0) then
        !bkgr1_test
           print *,'URBAN'
           !NAP
           beta(NA,1)     = beta(NA,1)   *0.75
           betamax(NA,1)  = betamax(NA,1)*0.75
           betamin(NA,1)  = betamin(NA,1)*0.75
           beta(NA,3)     = beta(NA,3)*0.35
           !AIT
           beta(AI,1)     = beta(AI,1)   *1.15
           betamin(AI,1)  = betamin(AI,1)*1.15
           betamax(AI,1)  = betamax(AI,1)*1.15
           beta(AI,3)     = beta(AI,3)*0.60
           !ACC
           beta(AS,1)     = beta(AS,1)   *1.5
           betamin(AS,1)  = betamin(AS,1)*1.5
           betamax(AS,1)  = betamax(AS,1)*1.5
           beta(AS,3)     = beta(AS,3)*0.37
           !COA
           beta(CS,3)     = beta(CS,3)*0.80

           ! change mode width
           betamin(AI,2)  = 1.90
           betamax(AI,2)  = 2.15
           !update mass min & max
           betamin(:,3)  = beta(:,3)*0.90

        else
        !xprs1,2_test
           print *,'EXHAUST'
           !NUC
           betamax(NU,1)  = 7.6 ! nm
           betamax(NU,2)  = 1.30
           betamin(NU,3)  = beta(NU,3)*0.45
           betamax(NU,3)  = beta(NU,3)*2.00
           !NAP
           beta(NA,1)     = beta(NA,1)   *0.95
           betamax(NA,1)  = betamax(NA,1)*0.95
           betamin(NA,1)  = betamin(NA,1)*0.95
           beta(NA,3)     = beta(NA,3)*0.68
           !AIT
           beta(AI,1)     = beta(AI,1)   *0.77
           betamin(AI,1)  = betamin(AI,1)*0.77
           betamax(AI,1)  = betamax(AI,1)*0.77
           beta(AI,3)     = beta(AI,3)*0.80
           !ACC
           beta(AS,1)     = beta(AS,1)   *1.05
           betamin(AS,1)  = betamin(AS,1)*1.05
           betamax(AS,1)  = betamax(AS,1)*1.05
           beta(AS,3)     = beta(AS,3)*0.75
           !COA
           beta(CS,3)     = beta(CS,3)*0.80

           ! change mode width
           betamin(AI,2)  = 1.85
           betamax(AI,2)  = 2.10
           betamin(AS,2)  = 1.80
           betamax(AS,2)  = 1.95
           betamin(CS,2)  = 1.85
           betamax(CS,2)  = 2.20
           !update mass min & max
           betamin(:,3)  = beta(:,3)*0.90

        endif

! 'PM10 set'
      else

        print *,'PM10 set'

        print *,'mstart ',mstart
        print *,'Mtot(AI) ',beta(AI,3)

        ! update MASS min & max
        betamin(:,3) = betamin(:,3)*0.08

        do m = mstart, nm
           betamax(m,3) = beta(m,3)*300.0
        enddo

        if(beta(AI,3)<600.0) then
        !arctic_test
        ! ("arctic example")
           print *,'MARINE'
           !NUC
           beta(NU,1)     = beta(NU,1)   *0.62
           betamax(NU,1)  = betamax(NU,1)*0.62
           betamin(NU,1)  = betamin(NU,1)*0.62
           beta(NU,3)     = beta(NU,3)*16.0
           !NAP
           beta(NA,1)     = beta(NA,1)   *0.87
           betamax(NA,1)  = betamax(NA,1)*0.87
           betamin(NA,1)  = betamin(NA,1)*0.87
           beta(NA,3)     = beta(NA,3)*3.0
           !AIT
           beta(AI,1)     = beta(AI,1)   *1.09
           betamin(AI,1)  = betamin(AI,1)*1.09
           betamax(AI,1)  = betamax(AI,1)*1.09
           beta(AI,3)     = beta(AI,3)*4.0
           !ACC
           beta(AS,1)     = beta(AS,1)   *0.9
           betamin(AS,1)  = betamin(AS,1)*0.9
           betamax(AS,1)  = betamax(AS,1)*0.9
           beta(AS,3)     = beta(AS,3)*6.5

           ! change mode width
           betamin(AI,2)  = 1.55
           betamax(AI,2)  = 1.75
           !update mass min & max
           betamin(:,3)  = beta(:,3)*0.20

        else if (beta(AI,3)<10000.0) then
        !
           print *,'URBAN'
           !nothing special



        else
        !xprs3_test (same as xprs2_test)
           print *,'EXHAUST'
           !NUC
           betamax(NU,1)  = 7.2 ! nm
           betamax(NU,2)  = 1.30
           beta(NU,3)     = beta(NU,3)*0.90
           betamin(NU,3)  = beta(NU,3)*0.40
           betamax(NU,3)  = beta(NU,3)*2.40
           !NAP
           beta(NA,1)     = beta(NA,1)   *0.76
           betamax(NA,1)  = betamax(NA,1)*0.76
           betamin(NA,1)  = betamin(NA,1)*0.76
           beta(NA,3)     = beta(NA,3)*0.63
           !AIT
           beta(AI,1)     = beta(AI,1)   *1.23
           betamin(AI,1)  = betamin(AI,1)*1.23
           betamax(AI,1)  = betamax(AI,1)*1.23
           beta(AI,3)     = beta(AI,3)*0.86
           !ACC
           beta(AS,1)     = beta(AS,1)   *1.05
           betamin(AS,1)  = betamin(AS,1)*1.05
           betamax(AS,1)  = betamax(AS,1)*1.05
           beta(AS,3)     = beta(AS,3)*1.45

           ! change mode width
           betamin(AI,2)  = 1.95
           betamax(AI,2)  = 2.10
           betamin(AS,2)  = 1.80
           betamax(AS,2)  = 1.95
           betamin(CS,2)  = 1.85
           betamax(CS,2)  = 2.20
           !update mass min & max
           betamin(:,3)  = beta(:,3)*0.90

        endif

      endif

      !print *,'mm2 min',betamin(:,3)
      !print *,'mm2 max',betamax(:,3)
      !print *,'mm2 est',beta(:,3)


!     Generic Nucleation Mode (mstart=2)
!     Add small mass in nucleation mode if zero
!      if (mstart==2) then
!         beta(NU,1) =  9.6  ! nm
!         beta(NU,2) =  1.30
!         beta(NU,3) =  0.45 * dataa(2,2)/3000.0
!      endif


! ** Check here the parameter overview
      print *,'gmd min',betamin(:,1)
      print *,'gmd max',betamax(:,1)
      print *,'gmd est',beta(:,1)

      print *,'sig min',betamin(:,2)
      print *,'sig max',betamax(:,2)
      print *,'sig est',beta(:,2)

      print *,'mm  min',betamin(:,3)
      print *,'mm  max',betamax(:,3)
      print *,'mm  est',beta(:,3)

   !   stop


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

         !do k=1,np1
         !  print *,'ini simplex',(p(m,k,j),j=1,np)
         !enddo
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

       if ( ((m==2).and.(kount.gt.maxcnt2))  &
             .or. ((m.ne.2).and.(kount.gt.maxcnt1))  ) then
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
            .or.((m==2).and.(kount.gt.maxcnt2))          &
            .or.((m.ne.2).and.(kount.gt.maxcnt1)) )      &
            then
      
        ! End of simplex

 56     print *,'  ===> error criterion satisfied'

        if (fe_log) then
          write(funit_log,54)(p(m,ilo,j),j=1,np)
          write(funit_log,*)
          write(funit_log,'(1X,A58)') 'FITDIS error criterion fulfilled'
          write(funit_log,*)
        endif

        print *,'Normal termination of simplex in mode ',m

        ! stop

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
        ! print particle number estimate of Simplex
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

! ***
! *** WRITE MODE DATA TO INAERO.DAT
! ***
        call output_fitaero(beta,dense,ntot,wamass,shrmode,  &
                            binlo,binup,m,mstart,            &
                            gmdout,sigout)  ! GMD and SIGMA as inaero.dat

        if(m.eq.nm) then
           print *,'-----------------------------------------------------------'
           do n=1,nm
             print *,'GMD(nm) SIGMA in',n,gmdout(n)*1.e9,sigout(n)
           enddo
           print *,'-----------------------------------------------------------'
        endif


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

         !print*,'hello start walking again'
         go to 31   ! start walking again

       endif

! *** End of simplex algorithm


 140   print *,'end mode',m


! MUST RECALCULATE NUMBER SIZE DISTRIBUTION
! WITH FINAL GMD AND SIGMA
! GMD AND SIGMA ARE ONLY KNOWN AFTER THE LAST MODE !!!

       if(m.eq.nm) then
        do n=mstart,nm
         do i= binlo(n),binup(n)

! *** 1_ Mass size distribution (dm*dlindp)
          !  massout = ( (beta(n,3)-wamass(n))                                &
            massout = ( beta(n,3)                                & 
                         /(SQRT(2.*pi)*dataa(i,1)*1.e-9*LOG10(sigout(n))))   &
                         *EXP(-0.5*(LOG(dataa(i,1)*1.e-9/gmdout(n))/LOG(sigout(n)))**2.)


            if (i.eq.binup(nm)) then
              massout = massout * (dataa(i,1)-dataa(i-1,1))*1.e-9
            else 
              massout = massout * (dataa(i+1,1)-dataa(i,1))*1.e-9
              if (n.eq.NA) massout = massout*0.7   ! correction for NAP 
            endif

! *** 2_ Volume concentration (dim.less) 
            volcout = massout  / (conv*dens(n))

! *** 3_ Number concentration
            numout(n,i)  = volcout / ( (pi/6.)*((dataa(i,1)*1.e-9))**3.0 ) ! #/m^3

! *** Normalize to mean mode diameter
            numout(n,i)  = numout(n,i) * (dataa(i,1)*1.e-9) / gmdout(n)

            numout(n,i)  = numout(n,i) * 1.e-6                             ! #/cm^3

            numout(n,i)  = max(numout(n,i),0.0) 

        !    if (i.eq.binup(nm)) then
        !      print *,'bin out',n,i,dataa(i,1),dataa(i,1)-dataa(i-1,1),massout,numout(n,i)
        !    else 
        !      print *,'bin out',n,i,dataa(i,1),dataa(i+1,1)-dataa(i,1),massout,numout(n,i)
        !    endif

! *** 4_ Write number to output file
            write (funit_out_sizedis, 1000) tab,tab,dataa(i,1),tab,tab,dataa(i,2),tab,tab,numout(n,i)

         enddo
        enddo
       endif


!     End of loop over aerosol modes (m)
      end do



!     Close all i/o files
      close(funit_in_sizedis)
      close(funit_out_inaero)
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
      deallocate(numout)
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
 1000 format(3(2A,F0.2))

      end program prog_fitaero
