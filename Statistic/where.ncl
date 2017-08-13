load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

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
	fsty	=	1979
	basyst	=	1979		;---set start base year---
	basyed	=	1998		;---set end base year---

	neof	=	8		;---set eof order---
	optEOF	=	True
	optETS	=	True
	optEOF@pcrit	=	0
	optEOF@jopt	=	0
	eps     =       1e-2		;---set min error---

	nlat	=	(latN-latS)/2
	nlon	=	(lonR-lonL)/2
	nyear	=	endy-strty+1

	name	=	"reconstructed_giss+np_"+latS+"_"+strty+"-"+endy+"_mode"+neof

;*********************
; add datas
;*********************

	dir	=	"./"
	f	=	addfile(dir+"annual_giss+np_30_1979-2012.nc","r")	;---add data---	
	p	=	f->air(:,:,(strty-fsty):(endy-fsty))
	rp	=	p(year|:,lat|:,lon|:)
printVarSummary(rp)
	lat	=	rp&lat
	lon	=	rp&lon

	num_lat	=	dimsizes(lat)
	num_lon	=	dimsizes(lon)

	ave	=	new((/num_lat,num_lon/),float)
	do i = 0,num_lat-1
	 do j = 0,num_lon-1

	  ave(i,j)  = avg(rp((1979-strty):(1998-strty),i,j))
	  rp(:,i,j) = rp(:,i,j)-ave(i,j)
	 end do
	end do


;*********************
; dineof interpolation
;*********************

	ave1	=	avg(rp)
	rp	=	rp-ave1

	wp	=	rp
	slp	=	wp({lat|latS:latN},{lon|lonL:lonR},year|:)  		;reorder data
printVarSummary(slp)
	slp0	=	slp
	slp	=	where(ismissing(slp),0.0,slp)


	do neval = 1,neof
	print("********************next begin*******************************")
	  slp1    =       slp
	  d       =       1.0
          n       =       0
	  eof          =       new((/neval,nlat,nlon/),float)
	  eof_ts       =       new((/neval,nyear/),float)
	  

	  do while (d.gt.eps)
 	   eof		=	eofunc(slp,neval,optEOF)
   	   eof_ts	=	eofunc_ts(slp,eof,optETS)
printVarSummary(eof)
printVarSummary(eof_ts)
	   
	   do nn = 0,neval-1
	   eof_ts(nn,:) = eof_ts(nn,:) + eof_ts@ts_mean(nn)
	   end do
	   slp	=	eof2data(eof,eof_ts)

	   slp	=	where(ismissing(slp0),slp,slp0)
 	   slp  =       where(ismissing(slp),0.0,slp)
	   slp1	=	where(ismissing(slp1),0.0,slp1)

	   d	=	max(abs(slp1-slp))
	   slp1	=	slp
	   n	=	n+1
	   print(n)
	   print(d)

	   if(n.gt.2000) then
	   break
	   end if

	  end do
	  delete(eof)
	  delete(eof_ts)
	end do

;*********************
; set basic year
;*********************

	num_lat	=	dimsizes(slp&lat)		
	num_lon	=	dimsizes(slp&lon)		

        do i = 0,num_lat-1
         do j = 0,num_lon-1
	   ave(i,j) = avg(slp(i,j,basyst-strty:basyed-strty))
          slp(i,j,:) = slp(i,j,:)-ave(i,j)
         end do
        end do

;*********************
; creat files
;*********************

	fc	=	addfile(name+"_temp.nc","c")
	fc->temp=	slp


end
