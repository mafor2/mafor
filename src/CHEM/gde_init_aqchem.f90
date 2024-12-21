! <gde_init_aqchem.f90 - A component of the Multicomponent
!                     Aerosol Dynamics Model MAFOR>
!*****************************************************************************! 
!* 
!*    Copyright (C) 2011-2023  Matthias Steffen Karl
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
!*****************************************************************************!
!*    All routines written by Matthias Karl
!* 
!*****************************************************************************!
module gde_init_aqchem

   use messy_mecca_kpp_Global

implicit none

public :: readaqphase



contains

  subroutine readaqphase(caq)
    !----------------------------------------------------------------------
    !     
    !  Read initial aqueous phase concentrations
    !  Aq. concentrations in coarse mode in molec/cm^3
    !
    !      author
    !      -------
    !      Dr. Matthias Karl
    !
    !      purpose
    !      -------
    !      Read initial aqueous phase concentrations
    !      Aqueous phase concentrations in coarse mode in molec/cm^3
    !
    !      interface
    !      ---------
    !      inaqchem.dat
    !
    !
    !      method
    !      ------
    !      read ascii file with space or tab separated entries
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

     real( dp),dimension(nspec), intent(in out) :: caq


! local
     INTEGER, PARAMETER                :: aqspec=14

     REAL( dp), dimension(aqspec)      :: ca0  

     integer               :: stat
     integer               :: gm
     integer               :: it,jk
! fix for tab reading
     character(len=180)    :: inchemline
     character(len=1)      :: tab = ACHAR(9)
     character(len=12)     :: dummyname
!----------------------------------------------------------------


        open(27,file='inaqchem.dat',status='old', iostat=stat)
! open error handling    
        if (stat.ne.0) then
           write(6,*) 'File inaqchem.dat cannot be opened !'
           stop
        end if
        do jk=1,aqspec
! fix for tab reading
           read(27,'(a)') inchemline
           do it=1,LEN(inchemline)
         ! replace each tab in the line with a space
             if( inchemline(it:it).EQ.tab) inchemline(it:it)=' '
           end do
           read(inchemline,*) dummyname, gm, ca0(jk)
! end fix for tab reading
        end do  


       close(27)


      ! initialize aqueous phase species concentrations
      ! ca0 of some aqueous phase species [mcl/cm3(air)]
      ! in the coarse (droplet) mode
      ! all in molec/cm^3(air)
      ! this works only for modal aqueous chemistry

! aqueous mode 1 (= AI mode)
        c(   ind_O2_a( 1)) = ca0( 1)
        c(  ind_OHm_a( 1)) = ca0( 2)
        c(   ind_Hp_a( 1)) = ca0( 3)
        c(   ind_O3_a( 1)) = ca0( 4)
        c(  ind_O2m_a( 1)) = ca0( 5)
        c(  ind_Clm_a( 1)) = ca0( 6)
        c(ind_Feppp_a( 1)) = ca0( 7)
        c(ind_HSO4m_a( 1)) = ca0( 8)                                                
        c(ind_SO4mm_a( 1)) = ca0( 9)
        c(ind_H2SO4_a( 1)) = ca0(10)
        c( ind_NH4p_a( 1)) = ca0(11)
        c(ind_HCO3m_a( 1)) = ca0(12)
        c( ind_NO3m_a( 1)) = ca0(13)  
        c(  ind_DOC_a( 1)) = ca0(14)

! aqueous mode 2 (= AS mode)
        c(   ind_O2_a( 2)) = ca0( 1)  
        c(  ind_OHm_a( 2)) = ca0( 2)       
        c(   ind_Hp_a( 2)) = ca0( 3) 
        c(   ind_O3_a( 2)) = ca0( 4)        
        c(  ind_O2m_a( 2)) = ca0( 5)
        c(  ind_Clm_a( 2)) = ca0( 6)        
        c(ind_Feppp_a( 2)) = ca0( 7)
        c(ind_HSO4m_a( 2)) = ca0( 8)                                                       
        c(ind_SO4mm_a( 2)) = ca0( 9)
        c(ind_H2SO4_a( 2)) = ca0(10) 
        c( ind_NH4p_a( 2)) = ca0(11)
        c(ind_HCO3m_a( 2)) = ca0(12) 
        c( ind_NO3m_a( 2)) = ca0(13)  
        c(  ind_DOC_a( 2)) = ca0(14)

! aqueous mode 3 (= CS mode)
        c(   ind_O2_a( 3)) = ca0( 1)  
        c(  ind_OHm_a( 3)) = ca0( 2)       
        c(   ind_Hp_a( 3)) = ca0( 3) 
        c(   ind_O3_a( 3)) = ca0( 4)        
        c(  ind_O2m_a( 3)) = ca0( 5)
        c(  ind_Clm_a( 3)) = ca0( 6)        
        c(ind_Feppp_a( 3)) = ca0( 7)
        c(ind_HSO4m_a( 3)) = ca0( 8)                                                       
        c(ind_SO4mm_a( 3)) = ca0( 9)
        c(ind_H2SO4_a( 3)) = ca0(10) 
        c( ind_NH4p_a( 3)) = ca0(11)
        c(ind_HCO3m_a( 3)) = ca0(12) 
        c( ind_NO3m_a( 3)) = ca0(13)  
        c(  ind_DOC_a( 3)) = ca0(14)



  end subroutine readaqphase

end module gde_init_aqchem
