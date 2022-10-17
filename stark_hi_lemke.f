      FUNCTION STARK_HI_LEMKE (XLAM, XLAMCENTER, T, ED,
     >      STKTA_DWS, NWS, STKTA_TS, NTS, 
     >      STKTA_ES,  NES, STKTA_PS)       
C*****************************************************************
C     Function to obtain stark profile for a tabulated transition
C         as a function of electron density, temperature and wavelength
C input:
C     ED      : electron density per cubic centimeter
C     T       : temperature in K 
C     XLAM    : wavelength for which the profile should be calculated
C     XLAMCENTER : line center wavelength
C
C To prepare the work of this function, the broadening table for the 
C    considered line must have been stored into STKTA_PS by calling 
C    SUBROUTINE READ_H_STARKDATA (see there for more comments)
C The data array STKTA_PS has indices::   
C     NWS            Number of "scaled frequency"  points tabulated 
C     NTS            Number of temperature points tabulated 
C     NES            Number of electron densitiy points tabulated 
C
C output: stark profile at xlam in arbitrary units (!!)
C         - must be re-normalized after!
C******************************************************************

      DIMENSION STKTA_DWS(NWS), STKTA_TS (NTS), STKTA_ES (NES)
      DIMENSION STKTA_PS (NWS,NTS,NES)

      ELOG = ALOG10(ED)                                    
      TLOG = ALOG10(T) 
      DLAM = XLAM - XLAMCENTER

C***  Attention, only for symmetric profiles !!!!!!!
      DLAM_LOCAL = ABS(DLAM)

C***  The tables are over a re-scaled frequency 
C***   called "delta alpha" in the paper Lemke M.: 1997, A&A Suppl. 122, 285 
      DLAM_LOCAL = DLAM_LOCAL / (1.25E-9 * (10**(ELOG*2.0/3.0)))

C***  Find T interval (INDX_T, INDX_T+1)
      TLOG = MAX(STKTA_TS(1),MIN(STKTA_TS(NTS),TLOG))
      DO I = 2, NTS
         INDX_T=I-1
         IF (TLOG .LT. STKTA_TS(I)) EXIT
      END DO
      T_WGT=STKTA_TS(INDX_T+1)-STKTA_TS(INDX_T)
      T_WGT =(TLOG-STKTA_TS(INDX_T))/T_WGT

C***  Find electron density interval (INDX_ED, INDX_ED+1)
      ELOG = MAX(STKTA_ES(1),MIN(STKTA_ES(NES),ELOG))
      DO I = 2, NES                                     
         INDX_ED = I-1                                           
         IF (ELOG .LT. STKTA_ES(I)) EXIT
      END DO
      ED_WGT = (ELOG-STKTA_ES(INDX_ED))/
     >                 (STKTA_ES(INDX_ED+1)-STKTA_ES(INDX_ED))

C***  Find wavelength index interval (ML, ML+1)
      DLAM_LOCAL = MAX(STKTA_DWS(1),MIN(STKTA_DWS(NWS),DLAM_LOCAL))
      DO I = 2, NWS
         ML = I-1
         IF (DLAM_LOCAL .LT. STKTA_DWS(I)) EXIT
      END DO
      MLP = ML + 1
      WL_WGT = (DLAM_LOCAL - STKTA_DWS(ML)) 
     >       / (STKTA_DWS(MLP) - STKTA_DWS(ML))

C---------------------------------------------
C***  First at wavelength index ML
C***  Interpolation in temperature
      ST00 = STKTA_PS(ML,INDX_T  ,INDX_ED)
      ST01 = STKTA_PS(ML,INDX_T+1,INDX_ED)
      ST0  = T_WGT * (ST01 - ST00) + ST00

      ST00 = STKTA_PS(ML,INDX_T  ,INDX_ED+1)
      ST01 = STKTA_PS(ML,INDX_T+1,INDX_ED+1)
      ST1  = T_WGT * (ST01 - ST00) + ST00

C***  Interpolation in density
      STARK0 = ED_WGT*(ST1-ST0)+ST0
C---------------------------------------------
C***  second at wavelength index ML+1 alias MLP
C***  Interpolation in temperature
      ST00 = STKTA_PS(MLP,INDX_T  ,INDX_ED)
      ST01 = STKTA_PS(MLP,INDX_T+1,INDX_ED)
      ST0  = T_WGT * (ST01 - ST00) + ST00

      ST00 = STKTA_PS(MLP,INDX_T  ,INDX_ED+1)
      ST01 = STKTA_PS(MLP,INDX_T+1,INDX_ED+1)
      ST1  = T_WGT * (ST01 - ST00) + ST00

C***  Interpolation in density
      STARK1 = ED_WGT*(ST1-ST0)+ST0
C---------------------------------------------

C***  Interpolation in wavelength
      STARK = WL_WGT * (STARK1 - STARK0) + STARK0

      STARK_HI_LEMKE = 10.**STARK

      RETURN
      END
