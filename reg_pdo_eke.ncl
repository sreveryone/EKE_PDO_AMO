;************************************************************************
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
;***************************************************************************
begin
; 要做的是year_start~year_end期间的pdo指数与eke之间的回归
; 因为读取的eke资料是每周一次的观测资料，大小是7*52*（year_end-year_start），而pdo指数是每月的月平均资料，大小是12*（year_end-year_start）
; 因此，需要将时间维数调整为一致才行，这样才能调用regress函数。 具体做法是:
; Jan=31, Feb=28, Mar=31, Apr=30, May=31, Jun=30, Jul=31, Aug=31, Sep=30, Oct=31, Nov=30, Dec=31
;    
    file_path="/mnt/d/station_for_keyan/SLA/EKE_OBS_from_sla.nc"
    f=addfile(file_path, "r")
    ;print(f)
    lat_s=-30
    lat_n=70
    lon_w=120
    lon_e=240

    ; lat: 240个格点
    lat_south=floattoint((lat_s+90)/0.25) ; 400
    lat_north=floattoint((lat_n+90)/0.25) ; 640

    lon_west=floattoint(lon_w/0.25)
    lon_east=floattoint(lon_e/0.25)

    eke_init=f->EKE(:,lat_south:lat_north,lon_west:lon_east); 分辨率0.25 x 
    time=f->TREF_MSLA; 开始时间：1993.01.06 间隔为 7 day
    ;print(time)
    time_length=sizeof(time)/8
    ;print(time_length)
;---------------------------------------------------------------------------------;
;
;                           将每周的EKE观测数据处理为月平均的EKE数据
;
;__________________________________________________________________________________；
    day=(time-1)*7+6
    ; 因为当 mod(day(i),365).eq.0的时候，应该也属于12月，如若不处理，
    ; 在下文就会被归类到month(0)~month(1)，但实际上应该是归类到month(11)~month(12)
    ; 因此，需要对mod(day(i),365).eq.0的情况进行处理：
    do i = 0, time_length-1
        if (mod(day(i), 365).eq.0) then
            ; write branch
            day(i)=364
        end if
    end do
    ; 这样处理，day(i)会被归类到month(11)~month(12)
    k=1
    number=0
    flag=0

    month=ispan(0, 12, 1)
    ;print(month)
    month(0)=0
    month(1)=31
    month(2)=59
    month(3)=90
    month(4)=120
    month(5)=151
    month(6)=181
    month(7)=212
    month(8)=243
    month(9)=273
    month(10)=304
    month(11)=334
    month(12)=365
    ; tempt(lat_south:lat_north,lon_west:lon_east)
    tempt=eke_init(1,:,:)
    tempt=0.
    ;print(tempt)
    ;-------------------------------------------------------------------------------------------------------------------;
    ;
    ;                   eke 数组存放的是 将eke 每周数据归类到相应的月份并求月平均的资料
    ;                     时间范围：1993.01 ~ 2012.02
    ;--------------------------------------------------------------------------------------------------------------------;

    eke=eke_init(0:229,:,:)

    eke!1 = "lat"
    eke&lat = eke_init&NBLATITUDES

    eke!2 = "lon"
    eke&lon = eke_init&NBLONGITUDES



    do while ( flag.lt.time_length); 当flag=time_length时 exit this loop
        do i = flag, time_length-1,1

            if (mod(k, 12).eq.0) then
                if ((mod(day(i),365).le.month(12)).and.(mod(day(i),365).gt.month(11))) then                   
                    tempt(:,:)=tempt(:,:)+eke_init(i,:,:)
                    number=number+1
                else
                    break ; 每个月对应的时间段内的日子寻找完毕，exit loop
                end if
            else
                if ((mod(day(i),365).le.month(mod(k, 12))).and.(mod(day(i),365).gt.month(mod(k-1,12)))) then                
                    tempt(:,:)=tempt(:,:)+eke_init(i,:,:)
                    number=number+1
                else
                    break; exit loop
                end if
            end if

        end do

        ;print(k)
        ;print(number)
        q=i-1
        ;print(q)
        ;print(day(q))


        eke(k-1,:,:)=tempt(:,:)/number
        tempt=0
        number=0
        k=k+1
        flag=i
        ;print(flag)
    end do
;---------------------------------------------------------------------------------------------;
;
;                               pdo指数
;
;---------------------------------------------------------------------------------------------;
    file_pdo="/mnt/d/station_for_keyan/"+"pdo/pdo_year_month_index.txt"
    fpdo=asciiread(file_pdo,-1,"string")

    ypdo=stringtointeger(str_get_field(fpdo,1,","))
    mpdo=stringtointeger(str_get_field(fpdo,2,","))
    pdo0=stringtofloat(str_get_field(fpdo,3,","))

    ptpdo=ind(ypdo.ge.1993.and.ypdo.le.2011) ; 筛选出1993~2011年的pdo
    
    pdo=pdo0(ptpdo)
    ;print(pdo)
;---------------------------------------------------------------------;
;                           标准化处理
;
;----------------------------------------------------------------------;
    ; eke中的缺省值处理
    missing_value_amo=-1e+34
    ;print(missing_value)
    do i = 0, 229
        do y = 0, lat_north-lat_south
            do x = 0, lon_east-lon_west
                ; write loop content
                if (ismissing(eke(i,y,x))) then
                    ; write branch
                    eke(i,y,x)=0.
                end if
            end do
        end do
    end do
    pdo_normalization=dim_standardize_n(pdo, 0, 0)
    eke_normalization=dim_standardize_n(eke, 0, 0)
    ;rc           = regCoef_n(pdo_normalization,eke_normalization(0:227,:,:),0,0)  ;_n choose the element
   rc           = regCoef_n(pdo,eke(0:227,:,:),0,0)  ;_n choose the element
   rc!0 = "lat"
   rc&lat = eke&lat
   rc!1 = "lon"
   rc&lon = eke&lon


  copy_VarCoords(eke_init(1,:,:), rc)                ; copy lat,lon coords
  rc@long_name="regression  coefficient"
  ;print(rc@tval)
  ;printVarSummary(rc)
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
  prob&lat = eke&lat   ; assign coordinate values to named dimensions
  prob&lon = eke&lon
  ;prob@long_name="probability"
  ;printVarSummary(prob)
;-------------------------------------------------------------------------------------------
;                                       出图
;--------------------------------------------------------------------------------------


    res                             = True
    res@gsnMaximize                 = True

    wks_type="png"
    wks=gsn_open_wks(wks_type,"/mnt/d/station_for_keyan/reg_pdo_eke_1993~2011_Pacific")
    wks_type@wkWidth = 2400
    wks_type@wkHeight =2400
    res@cnFillPalette  = "cmocean_balance";颜色设置

;---These are sample resources you might want to set
    res@cnFillOn          = True     ; 打开等值线填充
    res@cnLinesOn         = False   ; 打开等值线
    res@cnLineLabelsOn    = False    ; turn off line labels
;---如果用的是高分辨率的数据最好使用的一条绘图命令，可以加快绘图
    res@cnFillMode        = "RasterFill"
    res@trGridType        = "TriangularMesh"  
;---map range

  res@mpCenterLonF = 180
  res@mpMinLatF = lat_s
  res@mpMaxLatF = lat_n
  res@mpMinLonF = lon_w
  res@mpMaxLonF = lon_e
;---------------------------------------------------------------;
;                          
;                           等值线间隔 
;
;---------------------------------------------------------------;
res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF = -600
res@cnMaxLevelValF = 600
res@cnLevelSpacingF = 100  ;间隔

  res@tiMainString = "regression of pdo_index & EKE during 1993-2011"

;---------------------避免报错-------------------
    res@gsnAddCyclic          = False  
;---------------------------------------------------
  plot1 = gsn_csm_contour_map(wks,rc,res)
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