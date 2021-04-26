!> \file mo_wqm_routing.f90

!> \brief Instream nitrate routing and associated assimilation uptake and denitrification 


!> \details This module provides a new routing method based on the modified river network \n
!>          information which was generated by intersecting the real river network and \n
!>          mRM model generated river network.

!> \authors Xiangqian Zhou
!> \date April 2021

MODULE mo_wqm_routing
   USE mo_kind, ONLY: dp,i4
IMPLICIT NONE

PUBLIC :: conc_routing_process

CONTAINS

! ------------------------------------------------------------------

  !     NAME
  !         read_routing_data

  !     PURPOSE
  !>        \brief Read the modified river network information

  !>        \details  Read the modified river network information by which was generated \n 
  !>                  intersecting the real river network and mRM model generated river network.

  !     CALLING SEQUENCE
  !         

  !     INTENT(IN)
  !>        None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None


  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiangqian Zhou
  !>        \date April 2021
  !>          
! 
subroutine read_L11_ID_mask(iBasin, L11_ID_mask)
   use mo_mrm_tools, only: get_basin_info_mrm
   use mo_mrm_global_variables, only: dirMorpho,level11
   use mo_read_spatial_data, only: read_spatial_data_ascii, read_header_ascii

   implicit none
   ! local variables
   integer(i4), intent(in)                  :: iBasin      ! basin 
   integer(i4)                              :: nrows11, ncols11     
   integer(i4), dimension(:), allocatable   :: L11_ID_mask,L11_ID_temp
   logical,     dimension(:,:), allocatable :: mask 
   logical,     dimension(:), allocatable   :: mask11
   character(256) :: fname
   integer(i4) :: funit
   integer(i4) :: nCells 
   integer(i4), dimension(:,:), allocatable    :: data_i4_2d  ! LinkID but not used 

   call get_basin_info_mrm(iBasin, 11, nrows11, ncols11) 
   allocate( mask(nrows11, ncols11) )
  ! allocate necessary variables at Level11
   ! allocate(level11%nrows       (iBasin))
   ! allocate(level11%ncols       (iBasin))
   ! allocate(level11%xllcorner   (iBasin))
   ! allocate(level11%yllcorner   (iBasin))
   ! allocate(level11%cellsize    (iBasin))
   ! allocate(level11%nodata_value(iBasin))

   ! Header (to check consistency)
   fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl("L11_mask.asc"))
   funit = 11
   call read_header_ascii(trim(fName), funit,   &
        level11%nrows(iBasin),     level11%ncols(iBasin), level11%xllcorner(iBasin), &
        level11%yllcorner(iBasin), level11%cellsize(iBasin), level11%nodata_value(iBasin))


   ! fName = trim(adjustl(dirMorpho(iBasin))) // trim(adjustl("L11_mask.asc"))
       ! 
    call read_spatial_data_ascii(trim(fName),funit, &
            level11%nrows(iBasin),     level11%ncols(iBasin), level11%xllcorner(iBasin),&
            level11%yllcorner(iBasin), level11%cellsize(iBasin),data_i4_2d, mask)
    ! read linkID on L11 
     nCells = size(data_i4_2d, dim=1)*size(data_i4_2d, dim=2)
     allocate ( L11_ID_temp(nCells) )
     L11_ID_temp(:) = RESHAPE( data_i4_2d, (/nCells/) )
       
     !create mask on L11
     allocate ( mask11(nCells))
     mask11(:)= RESHAPE( mask, (/nCells/) )

     ! create linkID
     allocate ( L11_ID_mask(nCells) )
     L11_ID_mask(:) = pack(L11_ID_temp(:), mask11(:))

    ! free memory
       deallocate(data_i4_2d,mask)
   end subroutine read_L11_ID_mask
  
! ------------------------------------------------------------------

  !     NAME
  !         conc_routing_process

  !     PURPOSE
  !>        \brief water qulity processes within routing. 

  !>        \details calculates river water concentration including in-stream processes including mixture in each node and 
  !>        denitrification and assimilatory in each reach(link).         

  !     INTENT(IN)
  !>        \param[in] "integer(i4)                     :: nNodes"
  !>        \param[in] "integer(i4)                     :: nLinks"
  !>        \param[in] "integer(i4)                     :: TS"                simulation time step [h] 
  !>        \param[in] "integer(i4)                     :: nsubstances"
  !river networks (links)
  !>        \param[in] "integer(i4), dimension(:)       :: nPerm"             basin routing order
  !>        \param[in] "integer(i4), dimension(:)       :: nLink_from"        from node
  !>        \param[in] "integer(i4), dimension(:)       :: nLink_to"          to node
  !>        \param[in] "real(dp), dimension(:)          :: nLink_length"      length of each reach(link)
  !variables from hydrological model
  !>        \param[in] "real(dp), dimension(:)          :: nNode_qOUT"        total outflow from L11 cells
  !>        \param[in] "real(dp), dimension(:)          :: nNode_qTR"      outflow from each reach, only IT=2 was called by WQM
  !variables for in-stream WQ processes
  !>        \param[in] "real(dp), dimension(:,:)        :: nNode_concOUT"     conc. of total outflow from L11 cells  
  !>        \param[in] "real(dp), dimension(:)          :: rivertemp11"       reach temperature
  !>        \param[in] "real(dp), dimension(:)          :: pardeniratew"      denitrification rate in aquatic 
  !>        \param[in] "real(dp), dimension(:),         :: parautouptkrate"   potiential autotrophic uptake mgNm-2d-1
  !>        \param[in] "real(dp), dimension(:)          :: parprodrate"       assimilatory rate in aquatic   

  !     INTENT(INOUT)
  !>        \param[inout] "real(dp), dimension(:,:)        :: nNode_interload"    variables to store input load at each node
  !>        \param[inout] "real(dp), dimension(:,:)        :: nNode_concTIN"      conc. of total inflow of each nodes
  !>        \param[inout] "real(dp), dimension(:)          :: nLink_riverbox"     water volumn of each reach (link)
  !>        \param[inout] "real(dp), dimension(:,:)        :: nLink_criverbox"    conc. in the water volumn
  !>        \param[inout] "real(dp), dimension(:)          :: nLink_yravg_q"      yearly average discharge of each reach
  !results
  !>        \param[inout] "real(dp), dimension(:,:)   :: nNode_concMod"      conc. at each node (ouput variable!)  
  !>        \param[inout] "real(dp), dimension(:)     :: aquatic_denitri"    amount of denitrification in aquatic system (instream)
  !>        \param[inout] "real(dp), dimension(:)     :: aquatic_assimil"    amount of assimilatory(primary production)

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None
  !         

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Sep 2016
  subroutine conc_routing_process(iBasin,nNodes, nLinks, no_day, TS, nPerm, nLink_from, nLink_to, nLink_length,  &
    nNode_qOUT, nNode_qTR, nNode_concOUT, nNode_interload, nNode_concTIN, nLink_riverbox, nLink_criverbox, nLink_yravg_q, &
    rivertemp11, nNode_concMod, aquatic_denitri, aquatic_assimil,pardeniratew, parautouptkrate, parprodrate, &
    nor_gr,rz_coeff, f_light)

use mo_wqm_global_variables,     only: L11_rivert_avg10, L11_rivert_avg20  !, L11_flightglobal state variables for moving average

implicit none
integer(i4),               intent(in)    :: iBasin
integer(i4),                     intent(in)    :: nNodes
integer(i4),                     intent(in)    :: nLinks
integer(i4),                     intent(in)    :: no_day     !number of day accounter in current year
integer(i4),                     intent(in)    :: TS         !simulation time step [h] 

!river networks (links)
integer(i4), dimension(:),       intent(in)    :: nPerm      ! basin routing order
integer(i4), dimension(:),       intent(in)    :: nLink_from ! from node
integer(i4), dimension(:),       intent(in)    :: nLink_to   ! to node
real(dp), dimension(:),          intent(in)    :: nLink_length ! length of each reach(link)
!variables from hydrological model
real(dp), dimension(:),          intent(in)    :: nNode_qOUT ! total outflow from L11 cells
real(dp), dimension(:),          intent(in)    :: nNode_qTR  ! outflow from each reach, for dim2 only IT was called by WQM
!variables for in-stream WQ processes
real(dp), dimension(:,:),        intent(in)    :: nNode_concOUT     ! conc. of total outflow from L11 cells  
real(dp), dimension(:,:),        intent(inout) :: nNode_interload   ! variables to store input load at each node
real(dp), dimension(:,:),        intent(inout) :: nNode_concTIN     ! conc. of total inflow of each nodes
real(dp), dimension(:),          intent(inout) :: nLink_riverbox    ! water volumn of each reach (link)
real(dp), dimension(:,:),        intent(inout) :: nLink_criverbox   ! conc. in the water volumn
real(dp), dimension(:),          intent(inout) :: nLink_yravg_q     ! yearly average discharge of each reach
real(dp), dimension(:),          intent(in)    :: rivertemp11       ! reach temperature
!results
real(dp), dimension(:,:),        intent(inout) :: nNode_concMod     ! conc. at each node (ouput variable!)  
real(dp), dimension(:),          intent(inout) :: aquatic_denitri   ! amount of denitrification in aquatic system (instream)
real(dp), dimension(:),          intent(inout) :: aquatic_assimil   ! amount of assimilatory(primary production)
!WQ parameters
real(dp), dimension(:),          intent(in)    :: pardeniratew      !denitrification rate in aquatic 
real(dp), dimension(:),          intent(in)    :: parautouptkrate   !potiential autotrophic uptake mgNm-2d-1
real(dp), dimension(:),          intent(in)    :: parprodrate       !assimilatory rate in aquatic  
real(dp),                        intent(in)    :: nor_gr            !normalised global radiation (daily)  
real(dp), dimension(:),          intent(inout)    :: rz_coeff          !rz shading coefficient dim = number of reaches
real(dp), dimension(:),          intent(inout)    :: f_light          !rz shading coefficient dim = number of reaches

!local
real(dp), dimension(nNodes)  :: temp_qTIN! store total inflow of each node, represents "netNode_qTIN" in hydrological routing
integer(i4)             :: sec_TS, DT      ! seconds per timestep, number of steps in one day
integer(i4)             :: k,i,j
real(dp)                :: width, depth, benthic_area   ![m] reach(link) width, depth and beneath area
real(dp)                :: newbox
integer(i4)             :: iNode, tNode
logical                                   :: flag
integer(i4), dimension(:), allocatable :: L1_L11_Id_mask  ! added by zhouxi

temp_qTIN = 0.0_dp
nNode_interload = 0.0_dp
!constant parameter  
sec_TS = TS* 3600_i4    ! [s]
DT = 24_i4 / TS

!! added by zhouxi

call read_L11_ID_mask(iBasin, L1_L11_Id_mask) 
    !allocate(L1_L11_Id_mask(size(mask11)))
    !L1_L11_Id_mask = pack(L1_L11_Id(:), mask11 ) ! added by zhouxi
write(123,*) L1_L11_Id_mask



do k=1, nLinks  
  i = nPerm(k)
  iNode = nLink_from(i)
  tNode = nLink_to(i)
  ! yearly average channel discharge (for empirial channel width and depth calculation)
  nLink_yravg_q(i) = nLink_yravg_q(i) + (nNode_qTR(iNode) - nLink_yravg_q(i)) / (90.0_dp )  
  ! empirial equations for channel width and depth
  width = 5.40_dp *nLink_yravg_q(i) **(0.50_dp)  ! from M. Rode (2016), EST
  depth = 0.27_dp *nLink_yravg_q(i) **(0.39_dp)  ! from J.A. Moody(2002), Earth Surface Processes and Landforms
  benthic_area = width * nLink_length(i)

  !total input at "iNode"
  temp_qTIN(iNode) = temp_qTIN(iNode) + nNode_qOUT(iNode)
  !conc. of total infow of node(at iNode)
  nNode_concTIN(iNode,:) = (nNode_interload(iNode,:) + &
             nNode_qOUT(iNode) * nNode_concOUT(iNode,:) * sec_TS) / (temp_qTIN(iNode) * sec_TS)
  newbox = nLink_riverbox(i) + temp_qTIN(iNode) * sec_TS
  

  !update conc. and water volumn of reach i
  nLink_criverbox(i,:) = (nLink_criverbox(i,:) * nLink_riverbox(i) + nNode_interload(iNode,:) + &
             nNode_qOUT(iNode) * nNode_concOUT(iNode,:) * sec_TS) / newbox          
  nLink_riverbox(i) = newbox
  !instream denitrification and primary production
  !10- and 20-day moving mean temperature of river water
  L11_rivert_avg10(i) = L11_rivert_avg10(i) + (rivertemp11(i) - L11_rivert_avg10(i)) / (10.0_dp )
  L11_rivert_avg20(i) = L11_rivert_avg20(i) + (rivertemp11(i) - L11_rivert_avg20(i)) / (20.0_dp )

   !!! check if there exist a stream then instream processes could occur, otherwise instream processes could not occur
  loop2: do j=1,size(L1_L11_Id_mask)
  if (i /= L1_L11_Id_mask(j)) then
    flag = .FALSE.
    cycle loop2
  else
    flag = .true.
    exit loop2
  end if
  end do loop2

  if (flag) then
  call instream_nutrient_processes(TS,no_day,i,  rivertemp11(i), L11_rivert_avg10(i), &
     L11_rivert_avg20(i), benthic_area, depth, nLink_riverbox(i), nLink_criverbox(i,:),aquatic_denitri(i), aquatic_assimil(i), &
     pardeniratew(i),parautouptkrate(i), parprodrate(i), nor_gr, rz_coeff(i), f_light(i))
  
  !update water volumn of reach(link) 
  nLink_riverbox(i) = nLink_riverbox(i) - nNode_qTR(iNode) * sec_TS
  !calculate the load from each upstream reach and store as input of the to_node (tNode) 
  nNode_interload(tNode,:) = nNode_interload(tNode,:) + nLink_criverbox(i,:) * nNode_qTR(iNode) * sec_TS
  !accumulate flow to the 'to node' of the current link, as the upstream inflow of 'to node'
  temp_qTIN(tNode) = temp_qTIN(tNode) + nNode_qTR(iNode)
  else
  !update water volumn of reach(link) without instream processes
  nLink_riverbox(i) = nLink_riverbox(i) - nNode_qTR(iNode) * sec_TS
  !calculate the load from each upstream reach and store as input of the to_node (tNode) 
  nNode_interload(tNode,:) = nNode_interload(tNode,:) + nLink_criverbox(i,:) * nNode_qTR(iNode) * sec_TS
  !accumulate flow to the 'to node' of the current link, as the upstream inflow of 'to node'
  temp_qTIN(tNode) = temp_qTIN(tNode) + nNode_qTR(iNode)
  end if 
  

end do
!*************************************
!accumulate at the outlet of catchment
!*************************************
iNode = nLink_from(nPerm(nLinks))
tNode = nLink_to(nPerm(nLinks))
temp_qTIN(tNode) = temp_qTIN(tNode) + nNode_qOUT(tNode)
nNode_concTIN(tNode,:) = (nNode_interload(tNode,:) + nNode_qOUT(tNode)* nNode_concOUT(tNode,:) * sec_TS) / &
                        (temp_qTIN(tNode) * sec_TS)

!variable for final concentration output
do i=1, nNodes
nNode_concMod(i,:) = nNode_concTIN(i,:)  
end do

end subroutine conc_routing_process

! ------------------------------------------------------------------

  !     NAME
  !         instream_nutrient_processes

  !     PURPOSE
  !>        \brief In-stream processes in a reach(link). 

  !>        \details Calculates in-stream processes in a sepcific reach(link), including denitrification and assimilatory.
  !>                 

  !     INTENT(IN)
  !>        \param[in] "integer(i4)                    :: TS"            time step         
  !>        \param[in] "integer(i4)                    :: i"         from node id of this sepcific reach
  !>        \param[in] "real(dp)                       :: rivertemp11"   river water temperature
  !>        \param[in] "real(dp)                       :: rt_avg10"      10-day's average water temperature
  !>        \param[in] "real(dp)                       :: rt_avg20"      20-day'S average water temperature
  !>        \param[in] "real(dp)                       :: benthic_area"  [m2]reach benthic area
  !>        \param[in] "real(dp)                       :: depth"         [m]average reach depth 
  !>        \param[in] "real(dp)                       :: pardeniratew"  parameter
  !>        \param[in] "real(dp)                       :: parautouptkrate"  parameter
  !>        \param[in] "real(dp)                       :: parprodrate"   parameter
  !     INTENT(INOUT)
  !>        \param[inout] "real(dp)                       :: riverbox"    river water volumn in each reach  
  !>        \param[inout] "real(dp), dimension(:)         :: criverbox"   conc. in water wolumn 


  !     INTENT(OUT)
  !>        \param[out] "real(dp)                       :: aqdenitri"   amount of N denitrified in river water 
  !>        \param[out] "real(dp)                       :: aqassimil"   amount of N assimilated in river water

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None
  !         

  !     LITERATURE
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Sep 2016
  !>    Modified
  !>        X. Yang Jun 2017 enabled different time-step (hourly)

subroutine instream_nutrient_processes(TS,no_day,i, rivertemp11, rt_avg10, rt_avg20, benthic_area, &
           depth, riverbox, criverbox, aqdenitri, aqassimil, pardeniratew,parautouptkrate, parprodrate,nor_gr, rz_coeff, f_light)

use mo_wqm_shadingeffect,      only: rz_shading_coeff
use mo_wqm_global_variables,   only: GR_file_exist

implicit none

integer(i4),                    intent(in)     :: TS
integer(i4),                    intent(in)     :: no_day     ! for light calculate -- shading effect
integer(i4),                    intent(in)     :: i
real(dp),                       intent(in)     :: rivertemp11
real(dp),                       intent(in)     :: rt_avg10
real(dp),                       intent(in)     :: rt_avg20
real(dp),                       intent(in)     :: benthic_area  
real(dp),                       intent(in)     :: depth    
real(dp),                       intent(inout)  :: riverbox  
real(dp), dimension(:),         intent(inout)  :: criverbox    
real(dp),                       intent(out)  :: aqdenitri    
real(dp),                       intent(out)  :: aqassimil    
real(dp),                       intent(in)     :: pardeniratew
real(dp),                       intent(in)     :: parautouptkrate
real(dp),                       intent(in)     :: parprodrate  
real(dp),                       intent(in)     :: nor_gr   ! for shading effect
real(dp),                       intent(inout)  :: rz_coeff ! coefficient for each reach
real(dp),                       intent(inout)  :: f_light  ! 5 days' moving average of rz_coeff


!local
real(dp)                 :: INpool, ONpool, tp_ave,aqassimil0
real(dp)                 :: f_temp, f_conc, f_tp, f_temp1, f_temp2
real(dp)                 :: DT
integer(i4)              :: k
!HYPE parameter
real(dp), parameter      :: tpmean = 0.005_dp
! constant variables
real(dp), parameter      :: halfsatINwater = 1.5  ! mg/l
real(dp), parameter      :: maxdenitriwater = 0.5_dp !
real(dp), parameter      :: maxprodwater = 0.5_dp  ! 
!real(dp), parameter      :: maxdegradwater = 0.5_dp  
real(dp), parameter      :: halfsatIPwater = 0.05_dp  
real(dp), parameter      :: activedepth = 1.0_dp  

DT = 24.0_dp / TS

INpool = riverbox * criverbox(1) / 1000.0_dp   !kg
f_temp = tempfactor(rivertemp11)
f_conc = criverbox(1) / (criverbox(1) + halfsatINwater )

!denitrification amount
aqdenitri = pardeniratew * f_temp * f_conc * benthic_area / DT
aqdenitri = min(maxdenitriwater*INpool, aqdenitri) 

!update pool and conc.
INpool = INpool - aqdenitri
if (riverbox >0.0_dp) then
criverbox(1) = INpool / riverbox * 1000.0_dp   !mg/l
else
criverbox(1) = 0.0_dp
end if
!primary production and mineralisation (inverse processes)
if (GR_file_exist) then
!##############################################
!##new method considering light availability ##
!##############################################
call rz_shading_coeff(i, no_day, nor_gr, rz_coeff)
f_light = f_light + (rz_coeff-f_light) / 5.0_dp

aqassimil =0.0_dp
aqassimil0 = 0.0_dp
if (rivertemp11 >= 0.0_dp) then
   aqassimil0 =  parautouptkrate* f_light * activedepth * benthic_area * 10E-6 /DT   ![kg]
   ! here the parautouptkrate (potential N autotrophic uptake = 300 mg N m-2 d-1) 
   aqassimil0 = min(0.9_dp *INpool*activedepth, aqassimil0 )  
   !mineralization (respiration) assumed part of autotrophic assimilation
   aqassimil = parprodrate*aqassimil0
   !aqassimil = min(maxprodwater*INpool*activedepth, aqassimil )
end if
else
!##method from HYPE##    
tp_ave = tpmean   ! since phorsphose is not implemented at the monment..
f_tp = tp_ave / (tp_ave + halfsatIPwater)
ONpool = riverbox * criverbox(2) /1000.0_dp
if (rivertemp11 > 0.0_dp) then
   f_temp1 = rivertemp11 /20.0_dp
else
   f_temp1 = 0.0_dp
end if  
f_temp2 = (rt_avg10 - rt_avg20) / 5.0_dp  
f_temp = f_temp1 * f_temp2
aqassimil =0.0_dp
if (rivertemp11 >= 0.0_dp) then
   aqassimil = parprodrate * f_temp * f_tp * activedepth * depth * benthic_area/ DT !
   if (aqassimil > 0.0_dp) then
      aqassimil = min(maxprodwater*INpool*activedepth, aqassimil ) 
   else
      aqassimil = max(-maxprodwater*ONpool*activedepth, aqassimil) 
   end if
end if
end if
INpool = INpool - aqassimil
ONpool = ONpool + aqassimil
if (riverbox >0.0_dp) then
criverbox(1) = INpool / riverbox * 1000.0_dp   !mg/l
criverbox(2) = ONpool / riverbox * 1000.0_dp   !mg/l
else
criverbox(1) = 0.0_dp
criverbox(2) = 0.0_dp
end if 

end subroutine instream_nutrient_processes 

real(kind=dp) function tempfactor(c_temp)
  
  implicit none
  real(dp),               intent(in)  :: c_temp     ! temperature
  real(dp)                            :: f_temp     ! function value
  
  f_temp = 2**((c_temp - 20.0_dp) / 10.0_dp)
  if (c_temp < 5.0_dp) f_temp = f_temp * (c_temp / 5.0_dp)
  if (c_temp < 0.0_dp) f_temp = 0.0_dp
  tempfactor = f_temp
    
 
end function tempfactor

END MODULE mo_wqm_routing