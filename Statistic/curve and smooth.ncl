load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; This script is to plot 0(60)-90 time series and their 10 year smooth

begin

;*****************************************************
; set parameters
;*****************************************************
        latS    =       72.5		;---set start lat---
	ymax	=	1.5e6		;---set YMax Value---
	ymin	=	-1.5e6
	linethx	=	2.0		;---set line thickness---
        latN    =       90
        lonL    =       0
        lonR    =       360

        strty	=       1979
        endy    =       2015   
	fsty	=	1979
	bastrty	=	1979
	baendy	=	2015

	time	=	ispan(strty,endy,1)
	ntime	=	endy-strty+1

;*****************************************************
; add data
;*****************************************************
	plot    =       new(4,graphic)
        eps_name=       "arctic_ERA_0.75x0.75_temp_"+strty+"-"+endy

  	dir	=	"/home/hmj/cmip5/moisture/ERA/energyzonal/Latent/"
  	f	= 	addfile(dir+"lte_waf_ERA-Interim_0.75x0.75_mon_1979-2015.nc", "r")
	pz	=	month_to_annual(f->lte,1)
	po	=	pz((strty-fsty):(endy-fsty),:,:)
	printVarSummary(po)
	nlat	=	dimsizes(po&latitude)
	nlon	=	dimsizes(po&longitude)

;*****************************************************
; set anomalies
;*****************************************************
        ave     =       new((/nlat,nlon/),float)
        do i = 0, nlat-1
         do j = 0, nlon-1
          ave(i,j)	= avg(pz((bastrty-fsty):(baendy-fsty),i,j))
          po(:,i,j)	= po(:,i,j)-ave(i,j)
         end do
        end do
	p1	= po(year|:,{latitude|0:latN},{longitude|lonL:lonR})
	p2      = po(year|:,{latitude|latS:latN},{longitude|lonL:lonR})

;*****************************************************
; weigh_average
;*****************************************************
	lat1	=       p1&latitude
        lon1	=       p1&longitude
        wlat1	= 	cos(lat1*0.0174532)
        clat1	= 	tofloat(wlat1)
	avep1	=	wgt_areaave_Wrap(p1,clat1,1.0,1)			; p must be (...,lat,lon)
	printVarSummary(avep1)		; 0-90

	lat2    =       p2&latitude
        lon2    =       p2&longitude
        wlat2   =       cos(lat2*0.0174532)
        clat2   =       tofloat(wlat2)
        avep2   =       wgt_areaave_Wrap(p2,clat2,1.0,1)                        ; p must be (...,lat,lon)
        printVarSummary(avep2)		;67.5-90
	

;*****************************************************
; calculate 10years smooth
;*****************************************************
	sp1	= avep1
	sp2	= avep2
	do i = 1, ntime-10
	 sp1(i)	= avg(avep1(i:i+9))
	 sp2(i) = avg(avep2(i:i+9))
	end do
	do i = ntime-9, ntime-2
	 sp1(i)	= avg(avep1(i:))
	 sp2(i) = avg(avep2(i:))
	end do
	printVarSummary(sp1)
	printVarSummary(sp2)

;***************************************************** 
; create plot
;***************************************************** 
 	wks   = gsn_open_wks ("eps",eps_name)                 

 	res		= True                    
 	res@gsnDraw	= False                 ; don't draw
 	res@gsnFrame	= False                 ; don't advance frame
	res@vpHeightF	= 0.4			; change aspect ratio of plot
	res@vpWidthF	= 0.7
  	res@trYMaxF	= ymax
	res@trYMinF     = ymin
  	res@trXMinF	= strty              	; set x-axis minimum
	res@trXMaxF     = endy
	res@tiYAxisString	= "Latent Heat Anomaly ["+"J/kg]"      	; add an axis title 
	res@tiXAxisString       = "Year"
;	res@tiMainString 	= "Arctic Temperature Anomalies(>"+latS+"N)"
;	res@gsnLeftString       = strty+"-"+endy	
;	res@gsnLeftString	= "Corr= "+spc
	res@txFontHeightF	= 0.0195            			; change title font heights
	res@tiXAxisFontHeightF  = 0.016
        res@tiYAxisFontHeightF  = 0.016
	res@tmXBLabelFontHeightF =0.014                                 ;coordinate size
        res@tmYLLabelFontHeightF =0.014

;	res@tmXBMode   = "Explicit"
;	res@tmXBValues = ispan(1979,2015,2)
;	res@tmXBLabels = tostring(ispan(1979,2015,2))

 	res@xyLineThicknessF   	= linethx
        res@xyMonoLineColor     = False                                 ; want colored lines
        res@xyLineColors        = "Black"                      		; colors chosen
        res@xyDashPatterns      = 0.
	plot(0)	= gsn_csm_xy (wks,time,avep2,res)   	 

	res@xyLineColors        = "black"
        res@xyDashPatterns      = 1.
        plot(1) = gsn_csm_xy (wks,time,sp2,res)

;*****************************************************
; Draw
;*****************************************************
	overlay(plot(0),plot(1))
	draw(plot(0))
	draw(plot(1))
	frame(wks)
 
end
