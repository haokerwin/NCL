load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
;********************
;  set parameters
;********************
        latS    =       60		;---set plot year---
        latN    =       90
        lonL    =       0
        lonR    =       360

        strty   =       1900		;---set start year---
        endy    =       2005				
        fsty    =       1900		;---set first year---

	ntime	=	endy-strty+1
	time	=	ispan(strty,endy,1)
	eps_name	=	"Sig_num_"+latS+"_"+strty+"-"+endy+"_cmip5"	

;********************
;  set plot type
;********************
	wks     =       gsn_open_wks("eps",eps_name)
	gsn_define_colormap(wks,"BlueDarkRed18")
	plot1 	= 	new(3,graphic)
	plot2	=	new(3,graphic)

; resource for a certain plot
        res     =       True			
	res@gsnDraw     = False        ; don't draw yet
        res@gsnFrame    = False        ; don't advance frame yet
        res@cnFillOn    =       True
        res@cnLinesOn   =       False
        res@mpMinLatF   =       latS
;	res@gsnCenterString       = strty+"-"+endy
	res@gsnCenterStringFontHeightF	= .045
	res@gsnLeftStringFontHeightF  = .035
        res@gsnPolar    =       "NH"
	res@cnLineLabelsOn       = False
	res@lbLabelBarOn         = False
	res@gsnMaximize          = True
	res@cnLevelSelectionMode= "ExplicitLevels"
	res@cnLevels	= fspan(-.03, .03, 21)       	; 1900-2005
;	res@cnLevels    = fspan(-.20, .20, 17)		; 1998-2005
;	res@cnLevels	= fspan(-.04, .04, 21)		; 1900-1943
;	res@cnLevels    = fspan(-.06, .06, 21)		; 1943-1965

; resource for the total plot
	resP	=	True
	resP@gsnMaximize         = True         ; large format
        resP@gsnPanelLabelBar    = True         ; add common colorbar
;	resP@txString            = "Arctic Temperature Trend(>"+latS+"N)" 
	resP@gsnPolar            = "NH"
        resP@mpFillOn            = True
        resP@mpOutlineOn         = True
        resP@mpMinLatF     	 = latS
	resP@lbTitleString       = "~S~o~N~C/a"
        resP@lbTitlePosition      = "Right"                          ; title location
        resP@lbTitleDirection     = "Across"                          ; letter angle
        resP@lbTitleAngleF        = 0.                               ; title angle
        resP@lbTitleFontHeightF   = 0.01                              ; font height
        resP@lbLabelFontHeightF   = 0.010

; resource for the shade area
        opt     = True
        opt@gsnShadeFillType    = "pattern"             ; pattern fill
        opt@gsnShadeHigh        = 2                     ; use pattern #2
        opt@gsnShadeLow = 17                            ; use pattern #17

        opt@cnLevelSelectionMode= "ManualLevels"        ; set manual contour levels
        opt@cnMinLevelValF      = 0.00                  ; set min contour level
        opt@cnMaxLevelValF      = .08                   ; set max contour level
        opt@cnLevelSpacingF     = 0.001                 ; set contour spacing
        opt@cnInfoLabelOn       = False                 ; turn off info label
        opt@cnLinesOn           = False                 ; do not draw contour lines
        opt@cnLineLabelsOn      = False                 ; do not draw contour labels
        opt@cnFillScaleF        = 1.0                   ; add extra density

;********************
;  add file
;********************
	nume    = (/"a","b","c"/)
	dname	= (/"ES32","GISS1200km","Cowtan"/)
	dir1	= "/home/hmj/cmip5/reference_data/2x2data/2.3figure-spatial3analysisSig/"
	f1	= addfile(dir1+"tas_Amon_Ensemble32_historical_r1i1p1_185001-200512.nc", "r")
	po1	= f1->tas(strty-fsty:endy-fsty,:,:)
	p1	= po1(year|:,{LAT|latS:latN},{LON|lonL:lonR})				;annual ES32

	dir2	= "/home/hmj/cmip5/reference_data/2x2data/"
	fsty1	= 1880	;GISS
	fsty2	= 1850	; Cowtan
	endy1	= 2014
	name1	= "GISS1200km-ERSSTv4_2x2_bas80-00_mon_ori_"+fsty1+"-"+endy1+".nc"
	name2	= "Cowtan_2x2_bas80-00_mon_ori_"+fsty2+"-"+endy1+".nc"
	f2	= addfile(dir2+name1,"r")
	f3	= addfile(dir2+name2,"r")
	
        op2	= month_to_annual(f2->tas((strty-fsty1)*12:(endy-fsty1)*12+11,:,:),1)
        p2	= op2(year|:,{LAT|latS:latN},{LON|lonL:lonR})				;annual GISS
	op3	= month_to_annual(f3->tas((strty-fsty2)*12:(endy-fsty2)*12+11,:,:),1)
	p3	= op3(year|:,{LAT|latS:latN},{LON|lonL:lonR})				;annual Cowtan

;*****************************************************
; based on 8100
;*****************************************************
        num_lat = dimsizes(p1&LAT)
        num_lon = dimsizes(p1&LON)

        ave1	= new((/num_lat,num_lon/),float)
	ave2	= ave1
	ave3	= ave1
        do i = 0,num_lat-1
         do j = 0,num_lon-1
          ave1(i,j)	= avg(p1((1981-strty):(2000-strty),i,j))
          p1(:,i,j)	= p1(:,i,j)-ave1(i,j)
	
	  ave2(i,j)	= avg(p2((1981-strty):(2000-strty),i,j))
	  p2(:,i,j)	= p2(:,i,j)-ave2(i,j)

	  ave3(i,j)	= avg(p3((1981-strty):(2000-strty),i,j))
	  p3(:,i,j)     = p3(:,i,j)-ave3(i,j)
         end do
        end do

	printVarSummary(p1)
	printVarSummary(p2)
	printVarSummary(p3)

;*****************************************************
; assembly
;*****************************************************
	do ndata = 0,2,1	

	if(ndata.eq.0) then
	 slp	= p1
	end if

	if(ndata.eq.1) then
         slp    = p2
        end if

	if(ndata.eq.2) then
         slp    = p3
        end if

	slpp    = slp(year|:,{LAT|latS:latN},LON|:)      ;reorder
        regas   = regCoef(time,slpp(LAT|:,LON|:,year|:))

	tval    = onedtond(regas@tval , dimsizes(regas))
        df      = onedtond(regas@nptxy, dimsizes(regas)) - 2
        b       = tval    ; b must be same size as tval (and df)
        b       = 0.5
        prob    = betainc(df/(df+tval^2),df/2.0,b)

        copy_VarMeta(slpp(0,:,:),regas)
        copy_VarMeta(slpp(0,:,:),prob)

        regas   = smth9(regas,0.5,0.2,True)
        printVarSummary(regas)
        printVarSummary(prob)

	delete(res@cnLevels)
	res@cnLevels    = fspan(-.03, .03, 21)	
	res@gsnRightString      = ""
	res@gsnLeftString 	= nume(ndata)
	res@gsnCenterString	= dname(ndata)

        plot1(ndata) = gsn_csm_contour_map_polar(wks,regas,res)


	delete(res@cnLevels)
        res@lbLabelBarOn        = False
        res@cnInfoLabelOn       = False
        res@cnLineLabelsOn      = False
	res@gsnRightString      = ""
        res@gsnLeftString       = ""
        res@gsnCenterString     = ""
        delete(res@gsnCenterString)
        res@cnLevels            = fspan(-.10, .10, 21)
        plot2(ndata) = gsn_csm_contour(wks, prob, res )
        plot2(ndata) = gsn_contour_shade(plot2(ndata),0.01,500.,opt)    ; (plot,lowevl,highevl,opt)
        overlay(plot1(ndata),plot2(ndata))


	delete(slp)
	delete(slpp)
	delete(regas)
	delete(tval)
	delete(df)
	delete(b)
	delete(prob)

	end do


	gsn_panel(wks,plot1,(/1,3/),resP)



end




