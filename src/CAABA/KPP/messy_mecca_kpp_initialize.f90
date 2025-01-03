! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Initialization File
! 
! Generated by KPP-2.2.3_rs3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : messy_mecca_kpp_Initialize.f90
! Time                 : Sat Nov 30 15:26:45 2024
! Working directory    : /home/matthias/MAFOR_GIT/mafor_v218/src/CAABA/mecca
! Equation file        : messy_mecca_kpp.kpp
! Output root filename : messy_mecca_kpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE messy_mecca_kpp_Initialize

  USE messy_mecca_kpp_Parameters, ONLY: dp, NVAR, NFIX
  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Initialize - function to initialize concentrations
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Initialize ( )


  USE messy_mecca_kpp_Global

  USE messy_mecca_kpp_Parameters

  INTEGER :: i
  REAL(kind=dp) :: x

  CFACTOR = 1.000000e+00_dp

  x = (0.)*CFACTOR
  DO i = 1, NVAR
    VAR(i) = x
  END DO

  x = (0.)*CFACTOR
  DO i = 1, NFIX
    FIX(i) = x
  END DO

! constant rate coefficients
  RCONST(10) = 1.8e-12
  RCONST(12) = 1
  RCONST(23) = 3.5e-12
  RCONST(36) = 1.2e-14
  RCONST(37) = 1300
  RCONST(41) = 1.66e-12
  RCONST(50) = 1.2e-12
  RCONST(60) = 3e-14
  RCONST(73) = 1.2e-12
  RCONST(79) = 1.4e-10
  RCONST(82) = 5.2e-12
  RCONST(83) = 6e-14
  RCONST(85) = 3.6e-14
  RCONST(86) = 1e-10
  RCONST(87) = 1.7e-12
  RCONST(88) = 5e-12
  RCONST(89) = 5e-12
  RCONST(90) = 5e-12
  RCONST(91) = 1e-12
  RCONST(92) = 6e-11
  RCONST(98) = 1.3e-10
  RCONST(99) = 2.53e-14
  RCONST(100) = 2.5e-11
  RCONST(101) = 4.3e-11
  RCONST(111) = 7e-15
  RCONST(118) = 2.3e-12
  RCONST(131) = 4e-12
  RCONST(135) = 3e-14
  RCONST(299) = 7.6e-11
  RCONST(308) = 6.3e-15
  RCONST(336) = 4.29e-11
  RCONST(346) = 2.1e-11
  RCONST(347) = 2.14e-11
  RCONST(355) = 1.39e-11
  RCONST(356) = 1.73e-11
  RCONST(358) = 4.51e-12
  RCONST(359) = 2.49e-11
  RCONST(360) = 1.44e-10
  RCONST(372) = 9.82e-11
  RCONST(373) = 6.97e-11
  RCONST(375) = 1.36e-10
  RCONST(376) = 7.33e-11
  RCONST(415) = 3.2e-11
  RCONST(499) = 3.4e-13
  RCONST(524) = 3.5e-13
  RCONST(556) = 8.8e-11
  RCONST(559) = 8.8e-11
  RCONST(574) = 7.46e-11
  RCONST(578) = 8.33e-11
  RCONST(579) = 3.39e-11
  RCONST(586) = 1.63e-11
  RCONST(588) = 1.27e-11
  RCONST(591) = 1.4e-11
  RCONST(600) = 1.72e-12
  RCONST(601) = 4.8e-12
  RCONST(603) = 4.75e-13
  RCONST(606) = 1.2e-15
  RCONST(607) = 1e-14
  RCONST(608) = 1e-15
  RCONST(609) = 1.4e-12
  RCONST(610) = 1.22e-10
  RCONST(613) = 2.19e-11
  RCONST(614) = 3.68e-11
  RCONST(615) = 3.67e-11
  RCONST(623) = 3e-13
  RCONST(624) = 2.2e-19
  RCONST(625) = 4.45e-11
  RCONST(626) = 6.18e-12
  RCONST(627) = 4e-11
  RCONST(628) = 2.31e-11
  RCONST(636) = 5.68e-12
  RCONST(638) = 2e-18
  RCONST(639) = 5.2e-11
  RCONST(640) = 4.66e-11
  RCONST(641) = 3.7e-11
  RCONST(648) = 2.62e-11
  RCONST(649) = 2.45e-11
  RCONST(655) = 3.7e-11
  RCONST(657) = 4.32e-11
  RCONST(664) = 3.44e-11
  RCONST(665) = 1.16e-12
  RCONST(671) = 1.78e-11
  RCONST(677) = 2.29e-11
  RCONST(679) = 3.59e-12
  RCONST(683) = 5.53e-13
  RCONST(684) = 3e-12
  RCONST(685) = 6.78e-13
  RCONST(686) = 3e-12
  RCONST(720) = 1e-17
  RCONST(733) = 1e-17
  RCONST(740) = 1.03e-10
  RCONST(742) = 2.4e-17
  RCONST(747) = 2.65e-11
  RCONST(751) = 2.4e-17
  RCONST(769) = 2.52e-11
  RCONST(770) = 2.88e-11
  RCONST(772) = 2.52e-11
  RCONST(773) = 3.81e-11
  RCONST(780) = 9.7e-12
  RCONST(788) = 2.4e-17
  RCONST(818) = 2.8e-17
  RCONST(825) = 10000
  RCONST(839) = 4.829e-16
  RCONST(855) = 1.9e-11
  RCONST(866) = 2.52e-11
  RCONST(871) = 7.49e-11
  RCONST(872) = 6.65e-11
  RCONST(880) = 3.12e-13
  RCONST(884) = 1.01e-10
  RCONST(889) = 1.33e-10
  RCONST(891) = 9.23e-11
  RCONST(896) = 2.64e-11
  RCONST(902) = 1.1e-10
  RCONST(903) = 4.33e-11
  RCONST(904) = 7.55e-11
  RCONST(906) = 7.19e-11
  RCONST(920) = 3.79e-12
  RCONST(921) = 1.38e-11
  RCONST(928) = 4.26e-12
  RCONST(935) = 4.5e-12
  RCONST(936) = 1.27e-12
  RCONST(938) = 1.23e-12
  RCONST(941) = 1.72e-11
  RCONST(947) = 7.48e-11
  RCONST(948) = 5.43e-11
  RCONST(950) = 7.52e-11
  RCONST(957) = 1.5e-12
  RCONST(963) = 2.74e-11
  RCONST(972) = 5.44e-11
  RCONST(974) = 2e-18
  RCONST(975) = 6.2e-11
  RCONST(976) = 4.38e-11
  RCONST(977) = 1e-12
  RCONST(978) = 8e-19
  RCONST(979) = 6.9e-11
  RCONST(986) = 3.06e-11
  RCONST(988) = 7.09e-11
  RCONST(989) = 1.69e-11
  RCONST(990) = 1.21e-10
  RCONST(992) = 8.76e-13
  RCONST(993) = 4.44e-12
  RCONST(994) = 1e-14
  RCONST(1001) = 3.59e-12
  RCONST(1002) = 2.53e-11
  RCONST(1003) = 3.59e-12
  RCONST(1018) = 1.01e-11
  RCONST(1024) = 7.11e-12
  RCONST(1025) = 8.69e-11
  RCONST(1026) = 3.22e-12
  RCONST(1027) = 1.33e-11
  RCONST(1029) = 1.16e-10
  RCONST(1030) = 7.7e-11
  RCONST(1031) = 3.6e-11
  RCONST(1038) = 4.05e-11
  RCONST(1040) = 7.3e-11
  RCONST(1041) = 9e-14
  RCONST(1042) = 9e-13
  RCONST(1048) = 6.07e-11
  RCONST(1049) = 9.9e-11
  RCONST(1050) = 9.2e-18
  RCONST(1051) = 1e-10
  RCONST(1052) = 8.01e-11
  RCONST(1053) = 2.6e-12
  RCONST(1054) = 3.47e-12
  RCONST(1062) = 2e-18
  RCONST(1063) = 6.08e-11
  RCONST(1069) = 9e-13
  RCONST(1076) = 7.66e-11
  RCONST(1080) = 9.2e-11
  RCONST(1086) = 4.06e-11
  RCONST(1090) = 2.25e-15
  RCONST(1091) = 3e-14
  RCONST(1092) = 3.8e-12
  RCONST(1094) = 3e-13
  RCONST(1095) = 4.6e-12
  RCONST(1101) = 9.42e-11
  RCONST(1113) = 1.07e-10
  RCONST(1117) = 1.23e-10
  RCONST(1123) = 8.16e-11
  RCONST(1129) = 9.77e-11
  RCONST(1136) = 3.28e-11
  RCONST(1138) = 6.68e-11
  RCONST(1140) = 6.45e-11
  RCONST(1152) = 4.37e-11
  RCONST(1153) = 3.6e-12
  RCONST(1154) = 1.31e-10
  RCONST(1167) = 4.38e-11
  RCONST(1175) = 6.7e-11
  RCONST(1181) = 4.75e-12
  RCONST(1182) = 8.83e-13
  RCONST(1188) = 1.2e-10
  RCONST(1192) = 1.27e-11
  RCONST(1196) = 3.31e-11
  RCONST(1209) = 1.4e-11
  RCONST(1210) = 4.65e-11
  RCONST(1218) = 5e-18
  RCONST(1219) = 7.99e-11
  RCONST(1220) = 2.05e-11
  RCONST(1221) = 6.03e-12
  RCONST(1222) = 2.4e-15
  RCONST(1237) = 2.8e-17
  RCONST(1239) = 9.64e-11
  RCONST(1240) = 7.16e-11
  RCONST(1241) = 7.99e-11
  RCONST(1262) = 1.15e-10
  RCONST(1263) = 1.07e-10
  RCONST(1265) = 2.8e-12
  RCONST(1273) = 5.98e-11
  RCONST(1274) = 6.29e-11
  RCONST(1275) = 5.96e-11
  RCONST(1277) = 7.04e-11
  RCONST(1278) = 3.06e-11
  RCONST(1279) = 4.06e-11
  RCONST(1280) = 4.66e-12
  RCONST(1281) = 1.1e-12
  RCONST(1282) = 1.06e-12
  RCONST(1284) = 1e-12
  RCONST(1285) = 2.3e-11
  RCONST(1288) = 4.65e-11
  RCONST(1289) = 5.03e-12
  RCONST(1290) = 6.83e-12
  RCONST(1306) = 7.83e-15
  RCONST(1307) = 5.1e-14
  RCONST(1323) = 2.05e-10
  RCONST(1324) = 8.56e-11
  RCONST(1325) = 1.42e-10
  RCONST(1326) = 7.95e-11
  RCONST(1337) = 1.53e-12
  RCONST(1340) = 9.58e-11
  RCONST(1349) = 9.29e-11
  RCONST(1350) = 8.96e-11
  RCONST(1355) = 1.29e-11
  RCONST(1359) = 3.45e-11
  RCONST(1366) = 1.09e-11
  RCONST(1370) = 1.86e-11
  RCONST(1372) = 2.63e-11
  RCONST(1381) = 9.65e-12
  RCONST(1382) = 6.57e-12
  RCONST(1383) = 2.96e-12
  RCONST(1385) = 3.04e-12
  RCONST(1390) = 1.62e-11
  RCONST(1391) = 1.84e-12
  RCONST(1392) = 3.94e-12
  RCONST(1398) = 3.61e-11
  RCONST(1399) = 2.56e-11
  RCONST(1405) = 8.35e-11
  RCONST(1406) = 4.96e-11
  RCONST(1407) = 4.01e-12
  RCONST(1408) = 1.01e-12
  RCONST(1409) = 2.61e-12
  RCONST(1410) = 9.32e-12
  RCONST(1411) = 3.9e-16
  RCONST(1412) = 7e-12
  RCONST(1413) = 1.2e-16
  RCONST(1414) = 1.5e-12
  RCONST(1415) = 1.7e-17
  RCONST(1416) = 5.8e-11
  RCONST(1421) = 6.16e-11
  RCONST(1427) = 6.16e-11
  RCONST(1432) = 2.88e-12
  RCONST(1434) = 1.3e-11
  RCONST(1439) = 1.05e-11
  RCONST(1444) = 2.05e-11
  RCONST(1445) = 2.64e-11
  RCONST(1453) = 6.6e-12
  RCONST(1454) = 1.02e-11
  RCONST(1462) = 2.69e-11
  RCONST(1463) = 3e-11
  RCONST(1465) = 2.52e-11
  RCONST(1473) = 7.29e-12
  RCONST(1474) = 1.55e-11
  RCONST(1479) = 2.63e-11
  RCONST(1480) = 3.07e-12
  RCONST(1482) = 1.2e-15
  RCONST(1483) = 1e-14
  RCONST(1484) = 1e-15
  RCONST(1488) = 1.04e-11
  RCONST(1490) = 6.77e-12
  RCONST(1491) = 2.917e-11
  RCONST(1492) = 1.89e-12
  RCONST(1493) = 1.41e-12
  RCONST(1494) = 2.917e-11
  RCONST(1495) = 1.52e-15
  RCONST(1502) = 2.77e-11
  RCONST(1503) = 4.29e-12
  RCONST(1504) = 6.46e-11
  RCONST(1508) = 1e-11
  RCONST(1510) = 2e-14
  RCONST(1520) = 3.66e-12
  RCONST(1521) = 6.65e-12
  RCONST(1522) = 9.73e-12
  RCONST(1527) = 2.75e-11
  RCONST(1528) = 2.25e-11
  RCONST(1533) = 8.01e-11
  RCONST(1534) = 7.03e-11
  RCONST(1537) = 1.2e-15
  RCONST(1538) = 1e-14
  RCONST(1539) = 1e-15
  RCONST(1542) = 1.2e-15
  RCONST(1543) = 1e-14
  RCONST(1544) = 1e-15
  RCONST(1549) = 5.47e-11
  RCONST(1550) = 5.47e-11
  RCONST(1564) = 1.33e-11
  RCONST(1565) = 4.7e-12
  RCONST(1581) = 2.51e-12
  RCONST(1587) = 9.58e-12
  RCONST(1588) = 9.16e-13
  RCONST(1589) = 9.16e-13
  RCONST(1601) = 9.1e-12
  RCONST(1608) = 2.51e-12
  RCONST(1615) = 2.51e-12
  RCONST(1616) = 5.67e-11
  RCONST(1617) = 2.6e-15
  RCONST(1633) = 5.9e-11
  RCONST(1637) = 8e-11
  RCONST(1641) = 2.1e-12
  RCONST(1642) = 1.6e-10
  RCONST(1645) = 5e-11
  RCONST(1651) = 2.1e-10
  RCONST(1652) = 5e-12
  RCONST(1654) = 1e-10
  RCONST(1660) = 1.5e-12
  RCONST(1661) = 9e-12
  RCONST(1663) = 1.22e-12
  RCONST(1664) = 2e-12
  RCONST(1665) = 3.4e-17
  RCONST(1669) = 2.4e-15
  RCONST(1679) = 7e-14
  RCONST(1680) = 7e-14
  RCONST(1681) = 7e-14
  RCONST(1682) = 7e-14
  RCONST(1683) = 5e-13
  RCONST(1684) = 3e-14
  RCONST(1685) = 500000
  RCONST(1686) = 5e-12
  RCONST(1687) = 1e-13
  RCONST(1688) = 600000
  RCONST(1689) = 5e-14
  RCONST(1690) = 1e-11
  RCONST(1691) = 5.7e-12
  RCONST(1692) = 1.9e-11
  RCONST(1693) = 1.5e-12
  RCONST(1694) = 1.8e-13
  RCONST(1697) = 1.5e-14
  RCONST(1698) = 3e-18
  RCONST(1699) = 6.1e-11
  RCONST(1700) = 3.2e-13
  RCONST(1703) = 7.7e-18
  RCONST(1704) = 3e-12
  RCONST(1705) = 1.1e-10
  RCONST(1706) = 3.4e-12
  RCONST(1707) = 8.5e-13
  RCONST(1708) = 8.5e-13
  RCONST(1709) = 8e-13
  RCONST(1710) = 1.1e-11
  RCONST(1711) = 2e-17
  RCONST(1714) = 5.5e-12
  RCONST(1716) = 1e-13
  RCONST(1717) = 1e-15
  RCONST(1718) = 1e-15
  RCONST(1719) = 2.5e-13
  RCONST(1720) = 2.2e-11
  RCONST(1721) = 2.5e-13
  RCONST(1722) = 3e-15
  RCONST(1724) = 1e-11
  RCONST(1725) = 3e-12
  RCONST(1726) = 170
  RCONST(1727) = 5.5e-12
  RCONST(1728) = 1e-11
  RCONST(1729) = 2e-12
  RCONST(1730) = 170
  RCONST(1731) = 5.5e-12
  RCONST(1732) = 5.89e-12
  RCONST(1733) = 0.0115
  RCONST(1734) = 0.0115
  RCONST(1736) = 5e-12
  RCONST(1737) = 1.5e-12
  RCONST(1738) = 10
  RCONST(1739) = 2.24e-14
  RCONST(1740) = 7.03e-11
  RCONST(1741) = 1.11e-11
  RCONST(1743) = 1.4e-12
  RCONST(1751) = 4e-13
  RCONST(1752) = 1.2e-11
  RCONST(1753) = 3.3e-10
  RCONST(2096) = 7.9e-12
  RCONST(2097) = 1e-12
  RCONST(2098) = 3.6e+08
  RCONST(2099) = 4.4e+06
  RCONST(2105) = 2.4e-17
  RCONST(2123) = 1.2e-19
  RCONST(2126) = 3e-13
  RCONST(2132) = 9.54e-20
  RCONST(2134) = 4.5e-12
  RCONST(2135) = 8.5e-12
  RCONST(2136) = 3e-12
  RCONST(2145) = 8.5e-12
  RCONST(2147) = 400000
  RCONST(2149) = 8.5e-12
  RCONST(2158) = 9.54e-20
  RCONST(2160) = 1.74e-11
  RCONST(2161) = 8.5e-12
  RCONST(2162) = 1.61e-11
  RCONST(2169) = 3.18e-13
  RCONST(2170) = 1.7e-13
  RCONST(2179) = 8.5e-12
  RCONST(2181) = 8.5e-12
  RCONST(2183) = 400000
  RCONST(2184) = 2.4e-15
  RCONST(2185) = 1.4e-11
  RCONST(2186) = 8.5e-12
  RCONST(2188) = 1.1e-10
  RCONST(2204) = 4e-11
  RCONST(2205) = 4e-11
  RCONST(2206) = 9.6e-06
  RCONST(2207) = 9.6e-06
  RCONST(2208) = 4e-11
  RCONST(2209) = 4e-11
  RCONST(2210) = 9.6e-06
  RCONST(2211) = 9.6e-06
  RCONST(2212) = 2e-11
  RCONST(2213) = 2e-11
  RCONST(2214) = 0.0005
  RCONST(2486) = 0
  RCONST(2965) = 0
  RCONST(3444) = 0
! END constant rate coefficients

! INLINED initializations

  rtol(:) = 1E-2_dp ! relative tolerance
  atol(:) = 1E1_dp  ! absolute tolerance
  IF ((ind_OH  >0).AND.(ind_OH  <=NVAR)) atol(ind_OH)  = 1._dp
  IF ((ind_NO3 >0).AND.(ind_NO3 <=NVAR)) atol(ind_NO3) = 1._dp
  IF ((ind_Cl  >0).AND.(ind_Cl  <=NVAR)) atol(ind_Cl)  = 1._dp
  IF ((ind_Br  >0).AND.(ind_Br  <=NVAR)) atol(ind_Br)  = 1._dp
  IF ((ind_O1D >0).AND.(ind_O1D <=NVAR)) atol(ind_O1D) = 1._dp

! End INLINED initializations

      
END SUBROUTINE Initialize

! End of Initialize function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE messy_mecca_kpp_Initialize

