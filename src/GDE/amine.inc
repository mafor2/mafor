!************     MEA data sheet    ******************

  REAL( dp),    parameter :: xmAMIN     = 6.108E-2  ![kg/mol]
  REAL( dp),    parameter :: fcAMIN     = 0.393
  integer, parameter :: vhAMIN          = 2    

  ! Henry coefficient (mol/kg atm^-1)
  !
  REAL( dp),    parameter :: henry_AMIN = 6.20E06


  ! accomodation coefficient
  ! 
  REAL( dp),    parameter :: alpha_AMIN = 0.009 

  ! surface tension
  ! From Fit to MEA, MDEA and their mixtures 
  ! a=0.01741 b=138.9998
  !
  REAL( dp),    parameter :: a_AMIN = 0.01741   
  REAL( dp),    parameter :: b_AMIN = 138.998


!******************************************************
