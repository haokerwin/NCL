load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

; plot curves using different Y labels 
; @author hmj

begin
;*****************************************************
; set parameters
;*****************************************************
	strty   = 1979
        endy    = 2015
	bastrty = 1979
        baendy  = 2015
        fsty    = 1979
        latS    = 72
	latN	= 90
	lonL    = 0
        lonR    = 360
	nyear   = endy-strty+1

	pi      = 3.1415
        a       = 6.371e6
	Lv	= 2501.*1000.   ;J/g    notice the unit kg and g
	nseason	= "DJF"

	cname	= "Low_AMT"	;choose the case name	"High_SAT" and "Low_SAT"
	vname	= "SAT"		;choose the variable name
	uname	= "~S~o~N~C"	;choose the unit name
	name	= cname+"_"+nseason+"_temp_AMT"+latS

;*****************************************************
; add data
;*****************************************************
	dir2	= "/home/hmj/data/"
        ft2     = addfile(dir2+"ERA-Interim_0.75x0.75_mon_1979-2015_sfc_MSE.nc", "r")
	tt2	= month_to_season(short2flt(ft2->t2m),nseason)			; tt3(time,latitude,longitude)

	ta2	= tt2((strty-fsty):(endy-fsty),:,:)

	to2	= ta2(time|:,{latitude|latS:latN},{longitude|lonL:lonR})

        printVarSummary(to2)
	nlat	= dimsizes(to2&latitude)
	nlon	= dimsizes(to2&longitude)

	thigh	= new((/15,nlat,nlon/),float)	;an array to record the high amt year data
	tlow	= new((/15,nlat,nlon/),float)	;an array to record the low amt year data

;***************************************************** 
; read the order
;*****************************************************
	dir3	= "/home/hmj/cmip5/moisture/ERA/6-composite/"
	index_name	= "sortAMT79-15"
	nrow    = numAsciiRow(dir3+index_name+".txt")
	ncol    = numAsciiCol(dir3+index_name+".txt")	
	ft3	= asciiread(dir3+index_name+".txt",(/nrow,ncol/),"float")
	ord	= ft3(:,0)-strty
	ordl	= ord(0:14)			;the order of low
	ordh	= ord(nyear-15:nyear-1)		;the order of high
	
	do ii=0,15-1,1
	 thigh(ii,:,:)	= to2(toint(ordh(ii)),:,:)
	 tlow (ii,:,:)	= to2(toint(ordl(ii)),:,:)
	end do
	printVarSummary(thigh)				;thigh(time15,latitude,longitude)
	
;*****************************************************
; average along time and plot
;*****************************************************
	sthig	= dim_avg_n_Wrap(thigh,0)-273.15	
	stlow	= dim_avg_n_Wrap(tlow,0)-273.15		;sthig(latitude,longitude)
	stdif	= sthig - stlow
	copy_VarMeta(thigh(0,:,:),sthig)
	copy_VarMeta(tlow (0,:,:),stlow)
	copy_VarMeta(thigh(0,:,:),stdif)

;*****************************************************
; average of the whole area
;*****************************************************
	lat	= thigh&latitude
	lon	= thigh&longitude
	wlat	= cos(lat*0.0174532)
	clat	= tofloat(wlat)
	avhi	= wgt_areaave_Wrap(sthig,clat,1.0,1)
	avlw	= wgt_areaave_Wrap(stlow,clat,1.0,1)
	avdf	= avhi - avlw
	print("The high average value is "+avhi)
	print("The low average value is "+avlw)
	print("The different average value is "+avdf)	

;*****************************************************
; set resources and plot figures
;*****************************************************
	wks     = gsn_open_wks("eps",name)
	gsn_define_colormap(wks,"temp_19lev")
	res	= True
	res@cnFillOn	= True
	res@cnLinesOn	= False
	res@mpMinLatF	= latS
	res@gsnPolar	= "NH"
	
;	res@tiMainString	= cname
	res@gsnCenterString	= cname
	res@gsnRightString	= vname
	res@gsnLeftString	= uname 
	res@gsnCenterStringFontHeightF  = .030
        res@gsnLeftStringFontHeightF    = .015
        res@gsnRightStringFontHeightF   = .015

	res@cnLevelSelectionMode	= "ExplicitLevels"
	res@cnLevels	= fspan(-40.,0.,11)

	plot = gsn_csm_contour_map_polar(wks,stlow,res)		;sthig for the high plot, and the stlow for the low plot.

end

