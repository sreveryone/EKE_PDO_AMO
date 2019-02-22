;加载库函数
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin
    filepath1="/mnt/d/station_for_keyan/SST/HadISST_sst.nc";1870.01~
    ;filepath2="/mnt/d/station_for_keyan/SLA/EKE_OBS_from_sla.nc"
    f1=addfile(filepath1,"r")
    ;f2=addfile(filepath2,"r")
    ;print(f1)  ;查看nc文件属性
;---------------------------------------------------------------------;
;                    sst的趋势:2003.01~2012.12
    mo_start=1+12*(1993-1870)-1
    mo_end=12+12*(2011-1870)-1;数组从0开始所以要减一
;---------------------------------------------------------------------;
    ; y是三维：(time | 1785)x(lat | 180)x(lon | 360)
    ; 循环求每个格点的线性趋势beta
    ; trend为(lat | 180)x(lon | 360)
    ;lat_n=89
    ;lat_s=-89
    ;lon_w=-179.5
    ;lon_e=179.5

    lat_n=70
    lat_s=10
    lon_w=-179
    lon_e=179
    lat_south=89-lat_s
    lat_north=89-lat_n
    lon_west=lon_w+179
    lon_east=lon_e+180
    ; (lat | 19:159) x (lon | 0:358): 141x359
    sst0=f1->sst(mo_start:mo_end,lat_north:lat_south,lon_west:lon_east)
    printVarSummary(sst0)
    tempt=sst0
    do i=0,179,1
        tempt(:,:,i+180)=sst0(:,:,i)
    end do
    do i = 180, 359,1
        tempt(:,:,i-180)=sst0(:,:,i)
    end do


    time=f1->time(mo_start:mo_end)
    y=tempt
    ;y=sst0
    printVarSummary(tempt)
    x=time  
    lat1=y&latitude
    lon1=y&longitude


    dims_x=0;时间维数
    dims_y=0;时间维数
    rc  =   regCoef_n(x, y, dims_x, dims_y) ; rc(nlat,mlon)  
    rc=rc*10000
    rc!0             =   "lat"
    rc!1             =   "lon"
    rc&lat           =   lat1
    rc&lon           =   lon1
    printVarSummary(rc)
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
  prob&lat = lat1   ; assign coordinate values to named dimensions
  prob&lon = lon1
;---------------------------------------------------------------;
;
;                         绘图参数参数设置
;
;---------------------------------------------------------------;
    res                             = True
    res@gsnMaximize                 = True
    wks=gsn_open_wks("png","/mnt/d/station_for_keyan/sst_trend_1993-2011_Pacific")
    ;wks=gsn_open_wks("png","/mnt/d/station_for_keyan/sst_trend_1870~2017")
    ;gsn_define_colormap(wks,"cmp_b2r");2
    ;gsn_define_colormap(wks,"hotcolr_19lev");
    
    res@cnFillPalette  = "cmocean_balance";颜色设置

;---These are sample resources you might want to set
    res@cnFillOn          = True     ; 打开等值线填充
 
    res@cnLinesOn         = True    ; 打开等值线
    res@cnLineLabelsOn    = False    ; turn off line labels
;---如果用的是高分辨率的数据最好使用的一条绘图命令，可以加快绘图
;    res@cnFillMode        = "RasterFill"
;    res@trGridType        = "TriangularMesh"    
res@gsnLeftString  = "sst_trend_1993-2011"
;---------------------------------------------------------------;
;                          
;                           设置地图范围 
;
;---------------------------------------------------------------;
res@mpMinLatF = lat_s
res@mpMaxLatF = lat_n
res@mpMinLonF = 120
res@mpMaxLonF = 240
;res@mpMinLonF = lon_w
;res@mpMaxLonF = lon_e
res@mpCenterLonF = 180
;---------------------------------------------------------------;
;                          
;                           等值线间隔 
;
;---------------------------------------------------------------;
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -3.5
res@cnMaxLevelValF = 3.5
res@cnLevelSpacingF = 0.5  ;间隔
;---------------------避免报错-------------------
    res@gsnAddCyclic          = False  
;---------------------------------------------------
;---------------------------------------------------------------;
;                          绘图 
;---------------------------------------------------------------;
    plot1 = gsn_csm_contour_map(wks, rc, res)
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

  res2@cnLevels = (/0.1/)                       ;; 设置显著性水平
  res2@cnFillPatterns = (/17,-1/)
  res2@cnFillColors = (/1,1/)
  res2@gsnLeftString = " "
  plot2 = gsn_csm_contour(wks,prob,res2)
  overlay(plot1,plot2)

  draw(plot1)
  frame(wks)     
end