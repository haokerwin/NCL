load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;This routine tries to produce the annual ERA data from 1979-2012 as well 
;as ERA data which is masked at the same points as giss data from 1999-2012.

begin

;*********************
; set parameters
;*********************

	latS	=	30		;---set start lat---
	latN	=	90
	lonL	=	-180
	lonR	=	180

	strty	=	1979		;---set start year---
	endy	=	2012		;---set end year---
	fsty	=	1880

;*********************
; add datas
;*********************

	dir	=	"./"
	f1	=	addfile("2x2era_surface.nc","r")
	f	=	addfile("2x2Giss.air.2x2.250.mon.anom.comb.nc","r")		
	p	=	month_to_annual(f->air((strty-1880)*12:(endy-1880)*12+11,:,:),1)		;giss annual
	p1	=	month_to_annual(short2flt(f1->t2m((strty-1979)*12:(endy-1979)*12+11,:,:)),1) 	;ERA annual

	p1	=	lonFlip(p1)
	p	=	lonFlip(p)

	rp	=	p1
	rp(1999-strty:endy-strty,:,:) =	mask(rp(1999-strty:endy-strty,:,:),ismissing(p(1999-strty:endy-strty,:,:)),False)
	;mask 99-12ERA

;*********************
; creat files
;*********************

	fc1	=	addfile("annual_era.nc","c")
	fc1->temp=	p1
	fc2	=	addfile("masked_era.nc","c")
	fc2->temp=      rp
	fc3	=	addfile("annual_giss.nc","c")
	fc3->temp=      p


end
