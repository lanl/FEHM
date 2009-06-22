      subroutine stress_uncoupled(iflg) 
!**********************************************************************
!D1
!D1 PURPOSE
!D1
!D1 To manage the (fluid-stress) calculations
!D1 for for initial and final stress calculations
!D1
!**********************************************************************
!D2
!D2 REVISION HISTORY 
!D2
!D2 FEHM Version 2.20
!D2 
!D2 Initial implementation: Date 24-Oct-01, Programmer: George Zyvoloski
!D2
!D2 $Log:   /pvcs.config/fehm90/src/stressctr.f_a  $
!D2
!**********************************************************************
!D3 
!D3 REQUIREMENTS TRACEABILITY
!D3
!**********************************************************************
!D4 
!D4 SPECIAL COMMENTS AND REFERENCES
!D4 
!D4 Requirements from SDN: 10086-RD-2.20-00
!D4   SOFTWARE REQUIREMENTS DOCUMENT (RD) for the 
!D4   FEHM Application Version 2.20
!D4 
!**********************************************************************
c notes gaz
c 
c
      use comai
      use combi
      use comci
      use comdi
      use comdti
      use comei
      use comfi
      use comgi
      use comii
      use comki
      use comsi
      use comsplitts 
      use comxi
      use davidi
      
      implicit none

      integer iflg,i,ndummy,md,j,icoupl_save, ntty_save
      integer idof_stress_save,iad_save,itotal_save,itert_save
      character*10 macro1
      real*8, allocatable :: stressboun(:)
      integer, allocatable :: kq_dum(:)
      save icoupl_save, ntty_save, itotal_save, itert_save
      
      if(istrs.eq.0) return
      if(initcalc.eq.0) then 
         if(istrs_coupl.gt.0) return  
      endif
      
      
      if(iflg.eq.1) then
c     
c     initial stress calculations
c     using initial temperature and pressure
c     save displacements
c     
         if(istrs_coupl.eq.-2.or.istrs_coupl.eq.0.or.initcalc.ne.0) then
            idof_stress_save = idof_stress
            idof_stress = 3
            icoupl_save = istrs_coupl
            istrs_coupl = -99
            itotal_save = itotal
            iad_save = iad
            itert_save = itert
            call bnswer
            itert = itert_save 
            itotal = itotal_save
            iad_strs = iad
            iad = iad_save
            itotal_s = itotal_s + iad_strs
            idof_stress = idof_stress_save 
            istrs_coupl = icoupl_save
c     calculate total displacements
            call stressctr(10,0)
c     calculate stresses
            ntty_save = ntty 
            ntty = 2
c     calculate volume changes
            call stressctr(6,0)
c     update porosity
            call stressctr(-7,0)
            call stressctr(13,0)
         endif  
	 if(initcalc.ne.0) then
            initcalc = 0
c     leave porosity the same as inputted value	   
            ps = psini 
            vol_strain = 0.0
            vol_strain = 0.0
            if(icnl.ne.0) then
               du_ini = du
               dv_ini = dv
               duo = du
               dvo = dv
            else
               du_ini = du
               dv_ini = dv
               dw_ini = dw
               duo = du
               dvo = dv
               dwo = dw
            endif
	 endif
c     output diplacements and stresses
c     call stressctr(14,0)
 	 call stressctr(11,0)
	 ntty = ntty_save
         if(istrs_coupl.eq.0) then
c     just calculate stress and stop
            if (iout .ne. 0) then
               write(iout,*) 
               write(iout,*) 'stopping as requested in stress macro',
     &              ' (istrs_coupl.eq.0)'
            end if
            if(iptty.ne.0) then
               write(iptty,*)
               write(iptty,*)'stopping as requested in stress macro',
     &              ' (istrs_coupl.eq.0)'
            end if
         endif
      else if(iflg.eq.2) then
c     
c     final stress calculations
c     
         if(istrs_coupl.eq.-2.or.istrs_coupl.eq.-1) then
            itert_s = 0
            ntty_save = ntty
            ntty = 2
            icoupl_save = istrs_coupl
            istrs_coupl = -99
            itotal_save = itotal
            iad_save = iad
            itert_save = itert
            call bnswer
            itert = itert_save 
            itotal = itotal_save
            iad_strs = iad
            iad = iad_save
            itotal_s = itotal_s + iad_strs
            istrs_coupl = icoupl_save
            ntty = ntty_save
c     calculate total displacements
            call stressctr(10,0)
c     calculate stresses
            ntty_save = ntty
            ntty = 2
c     calculate volume changes
            call stressctr(6,0)
c     update porosity
            call stressctr(-7,0)
c     calculate stresses
            call stressctr(13,0)
c     output diplacements and stresses
c     call stressctr(14,0)
            call contr(1)
            call contr(-1)
            call stressctr(11,0)
            ntty = ntty_save
         endif          
      else if(iflg.eq.3) then
c     
c     stress calculations after each time step
c     
         if(istrs_coupl.eq.-3) then
            itert_s = 0
            ntty_save = ntty
            ntty = 2
            icoupl_save = istrs_coupl
            istrs_coupl = -100
            itotal_save = itotal
            iad_save = iad
            itert_save = itert
            call bnswer
            itert = itert_save 
            itotal = itotal_save
            iad_strs = iad
            iad = iad_save
            itotal_s = itotal_s + iad_strs
            istrs_coupl = icoupl_save
            ntty = ntty_save
c     calculate total displacements
c     call stressctr(10,0)
c     calculate stresses
            ntty_save = ntty
            ntty = 2
c     calculate volume changes
c     call stressctr(6,0)
c     calculate stresses
c     call stressctr(13,0)
c     reset order of equations
            call stress_combine(2,0)
            ntty = ntty_save 
         endif
      else if(iflg.ge.100) then
         
      endif

      return
      end
