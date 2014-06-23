	module comcomp
!***********************************************************************
! Copyright 2010 Los Alamos National Security, LLC  All rights reserved
! Unless otherwise indicated,  this information has been authored by an 
! employee or employees of the Los Alamos National Security, LLC (LANS),
! operator of the  Los  Alamos National  Laboratory  under Contract  No.
! DE-AC52-06NA25396  with  the U. S. Department  of  Energy.  The  U. S.
! Government   has   rights  to  use,  reproduce,  and  distribute  this
! information.  The  public may copy  and  use this  information without
! charge, provided that this  Notice and any statement of authorship are
! reproduced on all copies.  Neither  the  Government nor LANS makes any
! warranty,   express   or   implied,   or   assumes  any  liability  or
! responsibility for the use of this information.       
!***********************************************************************
!	phi	  Pressure in water-phase (comdi)
!	phil	  Pressure in oil-phase
!	phig	  Pressure in gas-phase
!	pcp	  Array for oil-gas capillary pressure (comdi)
!	pcpo	  Array for oil-water capillary pressure 
!	fl	  Saturation of oil (comco2)
!	fg	  Saturation of gas (comco2)
!	fw	  Saturation of water (comco2)
!	fk	  Saturation of kerogen
!	yo	  Mass fraction of oil in gas-phase
!	xo	  Mass fraction of oil in oil-phase
!	yg	  Mass fraction of gas in gas-phase
!	xg	  Mass fraction of gas in oil-phase
!	yw	  Mass fraction of water in gas-phase (comco2)
!	xw	  Mass fraction of water in oil-phase (comco2)
!	rl_l	  Array for relative permeability of oil-phase (comco2)
!	rl_g	  Array for relative permeability of gas-phase
!	rl_w	  Array for relative permeability of water-phase (comco2)
!	oil_prop	Array for oil properties
!	gas_prop	Array for gas properties
!	wat_prop	Array for water properties (comco2)
!	sto1	  Mass in accumulation term
!	sto1	  Energy in accumulation term
!	diw	  Advection in water-phase (comco2)
!	dil	  Advection in oil-phase (comco2)
!	div	  Advection in gas-phase (comci)
!	sk	  Water mass source/sink term (comdi)
!	qh	  Water energy source/sink term (comdi)
!	sko	  Oil mass source/sink term
!	qho	  Oil energy source/sink term
!	skg	  Gas mass source/sink term
!	qhg	  Gas energy source/sink term
!***********************************************************************

	real*8, allocatable :: pcrit(:)
	real*8, allocatable :: tcrit(:)
	real*8, allocatable :: rcrit(:)
	real*8, allocatable :: zcrit(:)
	real*8, allocatable :: ent0(:)
	real*8, allocatable :: hcap(:,:)
	real*8, allocatable :: accn(:)
	real*8, allocatable :: flcomp(:)
	real*8, allocatable :: mwt(:)
	real*8, allocatable :: bicof(:,:)

	real*8, allocatable :: phil(:)
	real*8, allocatable :: phig(:)
	real*8, allocatable :: phol(:)
	real*8, allocatable :: phog(:)
	real*8, allocatable :: pcpow(:)
	real*8, allocatable :: pcpog(:)
	real*8, allocatable :: fk(:)
	real*8, allocatable :: fok(:)
	real*8, allocatable :: rl_g(:)
	real*8, allocatable :: sko(:)
	real*8, allocatable :: qho(:)
	real*8, allocatable :: esko(:)
	real*8, allocatable :: skg(:)
	real*8, allocatable :: qhg(:)
	real*8, allocatable :: eskg(:)
	real*8, allocatable :: skk(:)
	real*8, allocatable :: qhk(:)
	real*8, allocatable :: oil_prop(:)
	real*8, allocatable :: gas_prop(:)
	
	real*8, allocatable :: diwg(:)
	real*8, allocatable :: diwk(:)
	real*8, allocatable :: divg(:)
	real*8, allocatable :: divk(:)
	real*8, allocatable :: dilg(:)
	real*8, allocatable :: dilk(:)
	
	real*8, allocatable :: dmgf(:)
	real*8, allocatable :: dmkf(:)
	real*8, allocatable :: degf(:)
	real*8, allocatable :: dekf(:)
	
	real*8, allocatable ::  dqg(:)
	real*8, allocatable ::  dqk(:)
	real*8, allocatable ::  deqg(:)
	real*8, allocatable ::  deqk(:)
	real*8, allocatable ::  dqhg(:)
	real*8, allocatable ::  dqhk(:)

	real*8, allocatable ::  denoili(:)
	real*8, allocatable ::  denoilh(:)
	real*8, allocatable ::  deneoili(:)
	real*8, allocatable ::  deneoilh(:)

	real*8, allocatable ::  dengasi(:)
	real*8, allocatable ::  dengash(:)
	real*8, allocatable ::  denegasi(:)
	real*8, allocatable ::  denegash(:)

	real*8, allocatable ::  denkeri(:)
	real*8, allocatable ::  denkerh(:)
	real*8, allocatable ::  denekeri(:)
	real*8, allocatable ::  denekerh(:)

	real*8, allocatable ::  pflowoil(:)
	real*8, allocatable ::  eflowoil(:)
	real*8, allocatable ::  pflowgas(:)
	real*8, allocatable ::  eflowgas(:)

	real*8, allocatable :: xoi(:)
	real*8, allocatable :: xg(:)
!	real*8, allocatable :: xw(:)
	real*8, allocatable :: yoi(:)
	real*8, allocatable :: yg(:)
!	real*8, allocatable :: yw(:)


	real*8, allocatable :: xoo(:)
	real*8, allocatable :: xog(:)
!	real*8, allocatable :: xow(:)
	real*8, allocatable :: yoo(:)
	real*8, allocatable :: yog(:)
!	real*8, allocatable :: yow(:)

	real*8, allocatable :: dentmp(:)
	real*8, allocatable :: mutmp(:)
	real*8, allocatable :: entmp(:)
	real*8, allocatable :: eosa(:)
	real*8, allocatable :: eosb(:)
	real*8, allocatable :: ddenp(:)
	real*8, allocatable :: ddent(:)
	real*8, allocatable :: dmudt(:)
	real*8, allocatable :: dmudp(:)
	real*8, allocatable :: dentt(:)
	real*8, allocatable :: dentp(:)

	real*8, allocatable :: cpr_zero(:,:)
	real*8, allocatable :: cpr_one(:,:)

	real*8, allocatable :: tkval(:)
	real*8, allocatable :: pkval(:)
	real*8, allocatable :: kval(:)
	integer, allocatable :: itind(:)
	integer, allocatable :: ipind(:)

	integer, allocatable :: kao(:)
	integer, allocatable :: kag(:)
	character*10, allocatable :: compin(:)

	integer ncomp, ioiltype, ioil, idof_oil, nkvt, nkvp, nkv

	end module comcomp
