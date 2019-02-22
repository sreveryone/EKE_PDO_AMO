;加载库函数
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin
    ;filepath1="/mnt/d/station_for_keyan/SST/HadISST_sst.nc";1870.01~
    filepath2="/mnt/d/station_for_keyan/SLA/EKE_OBS_from_sla.nc"
    ;f1=addfile(filepath1,"r")
    f2=addfile(filepath2,"r")
    ;print(f2)  ;查看nc文件属性  0.25 x 0.25
    ;start:  1993.01.06  interval = 7days
    ;---------------------------------------------------------------------;
    ;                    eke的趋势:day_start:day_end
    day_start= 52*(1993-1993)
    day_end=   52*(2011-1993);7天一周，一年52周
    ;---------------------------------------------------------------------;
    ; y   ：       (time | 1~997)x(lat | 721: -90~90)x(lon | 1440 0~359.75)
    ; 循环求每个格点的线性趋势beta
    ; trend为(lat | 721)x(lon | 1440)
    lat_n=70
    lat_s=10
    lon_w=120
    lon_e=240
    lat_north=1+tointeger((lat_n+90)/0.25) ; 
    lat_south=1+tointeger((lat_s+90)/0.25)   ; 
    lon_west=1+tointeger(lon_w/0.25)         ; 
    lon_east=1+tointeger(lon_e/0.25)
    ; (lat | 121:601) x (lon | 481:961): 481x481
    eke_init=f2->EKE(day_start:day_end,lat_south:lat_north,lon_west:lon_east)
    ;printVarSummary(eke_init)
    time=f2->TREF_MSLA(day_start:day_end)
    y=eke_init
    x=time  
    lat=eke_init&NBLATITUDES
    lon=eke_init&NBLONGITUDES

    dims_x=0;时间维数
    dims_y=0;时间维数
    rc  =   regCoef_n(x, y, dims_x, dims_y) ; rc(nlat,mlon)  
    ;rc2=rc2*10000
    rc!0             =   "lat"
    rc!1             =   "lon"
    rc&lat           =   lat
    rc&lon           =   lon
;-------------------------------------------------------
;           t检验回归rc是否显著
;-------------------------------------------------------
  tval = onedtond(rc@tval , dimsizes(rc))   ;t-statistic of rc
  df   = onedtond(rc@nptxy, dimsizes(rc)) - 2  ;自由度
  ;b = tval    ; b must be same size as tval (and df)
  ;b = 0.5
  prob  = student_t(tval, df)
  ;prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)
  prob!0   = "lat"    ; name dimensions
  prob!1   = "lon"
  prob&lat = lat   ; assign coordinate values to named dimensions
  prob&lon = lon
;---------------------------------------------------------------;
;
;                         绘图参数参数设置
;
;---------------------------------------------------------------;
    res                             = True
    res@gsnMaximize                 = True

        wks_type="png"
    wks=gsn_open_wks(wks_type,"/mnt/d/station_for_keyan/eke_trend_1993~2011_Pacific")
    wks_type@wkWidth = 2400
    wks_type@wkHeight =2400

    ;wks=gsn_open_wks("png","/mnt/d/station_for_keyan/sst_trend_1870~2017")
    ;gsn_define_colormap(wks,"cmp_b2r");2
    ;gsn_define_colormap(wks,"hotcolr_19lev");
    res@cnFillPalette  = "cmocean_balance";颜色设置

;---These are sample resources you might want to set
    res@cnFillOn          = True     ; 打开等值线填充
    res@cnLinesOn         = False   ; 打开等值线
    res@cnLineLabelsOn    = False    ; turn off line labels
;---如果用的是高分辨率的数据最好使用的一条绘图命令，可以加快绘图
    res@cnFillMode        = "RasterFill"
    res@trGridType        = "TriangularMesh"  
;--------------------------- --------------------------------------；
;                              加标题                              ;
;--------------------------------------------------------------------;
res@gsnLeftString  = "eke_trend_1993-2011"
;---------------------------------------------------------------;
;                          
;                           设置地图范围 
;
;---------------------------------------------------------------;
res@mpMinLatF = lat_s
res@mpMaxLatF = lat_n
res@mpMinLonF = lon_w
res@mpMaxLonF = lon_e
res@mpCenterLonF = 180
;---------------------------------------------------------------;
;                          
;                           等值线间隔 
;
;---------------------------------------------------------------;
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -1
res@cnMaxLevelValF = 1
res@cnLevelSpacingF = 0.1  ;间隔

;---------------------------------------------------------------;
;                          绘图 
;---------------------------------------------------------------;
    plot1 = gsn_csm_contour_map(wks, rc, res)
;---------------------避免报错-------------------
    res@gsnAddCyclic          = False  
;---------------------------------------------------
;------------------------------------------------------------------------;
;                               绘制显著性区域
;
;------------------------------------------------------------------------------;
  res2 = True
  res2@gsnDraw = False
  res2@gsnFrame = False
  res2@cnFillOn = True
  res@cnLinesOn = False
  res2@cnLineLabelsOn = False
  res2@cnInfoLabelOn = False
  res2@lbLabelBarOn = False
  res2@cnMonoFillPattern = False ; 测试不可缺少
  res2@cnLevelSelectionMode = "ExplicitLevels"

  res2@cnLevels = (/0.05/)                       ;; 设置显著性水平
  res2@cnFillPatterns = (/17,-1/)
  res2@cnFillColors = (/1,1/)
  ;res2@gsnLeftString = " "
  plot2 = gsn_csm_contour(wks,prob,res2)
  overlay(plot1,plot2)

  draw(plot1)
  frame(wks)     
end