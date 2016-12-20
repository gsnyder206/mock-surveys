;#### 
;#### repackage_mini_snapshot.pro
;#### This code written by Paul Torrey (ptorrey@cfa.harvard.edu)
;#### Writes a new simulation snapshot with the desired particles
;#### This was copied from Paul's cosmo-idl repository as a static optional version for the mock UDF code.
;#### (https://bitbucket.org/ptorrey/cosmo-idl)
;#### 

pro repackage_mini_snapshot,	x0,y0,z0,vx0,vy0,vz0,m0,id0,$
				x1,y1,z1,vx1,vy1,vz1,m1,id1,$
				x2,y2,z2,vx2,vy2,vz2,m2,id2,$
				x3,y3,z3,vx3,vy3,vz3,m3,id3,$
				x4,y4,z4,vx4,vy4,vz4,m4,id4,$
				x5,y5,z5,vx5,vy5,vz5,m5,id5,$
				icname=icname,$
				time_in=time_in,$
				redshift_in=redshift_in,$
				internal_energy_gas=internal_energy_gas,$
				rho_gas = rho_gas,$
				filename=filename,$
				halo_position = halo_position,$
				true_snapshot_format = true_snapshot_format,$
				met0=met0,met4=met4,$
				stellarage=stellarage,$
				snapshot_name = snapshot_name,$
                                ea=ea,$
                                nha=nha,$
                                smoothinglength=smoothinglength,$
                                sfr=sfr,$
				volume=volume,$
				sunrise_only=sunrise_only,$
				boxsize=boxsize, $
                                initmass=gfm_initial_mass

print,min(id0),max(id0)
;print,min(id1),max(id1)
print,min(id4),max(id4)


;help,x0,x1,x2,x3,x4,x5

x0 = reform(x0)
y0 = reform(y0)
z0 = reform(z0)
vx0 = reform(vx0)
vy0 = reform(vy0)
vz0 = reform(vz0)
m0 = reform(m0)
id0 = reform(id0)

if keyword_set(x1) then begin
x1 = reform(x1)
y1 = reform(y1)
z1 = reform(z1)
vx1 = reform(vx1)
vy1 = reform(vy1)
vz1 = reform(vz1)
m1 = reform(m1)
id1 = reform(id1)
endif


x4 = reform(x4)
y4 = reform(y4)
z4 = reform(z4)
vx4 = reform(vx4)
vy4 = reform(vy4)
vz4 = reform(vz4)
m4 = reform(m4)
id4 = reform(id4)


nbackground_base = 16l
nbackground = nbackground_base^3


if not keyword_set(snapshot_name) then snapshot_name = 'snapshot_000.hdf5'
if not keyword_set(filename) then filename = 'dump.dat'
if not keyword_set(halo_position) then halo_position=[0,0,0]
fout = filename

npart=lonarr(6)	
massarr=dblarr(6)*0.0
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L
npartall=lonarr(6)	
flag_cooling= 0L
num_files= 1L  
;;;BoxSize = 20000.0D
bytesleft=120
la=intarr(bytesleft/2)

if keyword_set(time_in) 	then time 	+= time_in
if keyword_set(redshift_in) 	then redshift 	+= redshift_in


;---- before dealing with particle numbers:
;--  1) Add background everywhere
;--  2) Remove background within actual data
;--  3) recount number of background cells used
;--  4) proceed as usual
;-------------------------------------------

N0 = n_elements(x0)	; This number will change once bg cells are added
bg_number = 0
;;;sunrise_only=1

if not keyword_set(sunrise_only) then begin


  local_pos  = fltarr(3,nbackground_base^3)
  local_vel  = fltarr(3,nbackground_base^3)
  local_rho  = fltarr(nbackground_base^3)
  local_u    = fltarr(nbackground_base^3)
  local_mass = fltarr(nbackground_base^3)

    for iii=0,nbackground_base-1 do begin
      x_pos = (0.5001+1.0*iii)*boxsize/(nbackground_base*1.0) - boxsize*0.5
      for jjj=0,nbackground_base-1 do begin
        y_pos = (0.5001+1.0*jjj)*boxsize/(nbackground_base*1.0) - boxsize*0.5
        for kkk=0,nbackground_base-1 do begin
          z_pos = (0.5001+1.0*kkk)*boxsize/(nbackground_base*1.0) - boxsize*0.5

 	  relative_halo_position = sqrt( x_pos^2 + y_pos^2 + z_pos^2 )

	  if relative_halo_position gt 3000.0 then begin	  
;;= kkk + jjj*nbackground_base + iii * nbackground_base^2

	    local_pos[0,bg_number] = x_pos
	    local_pos[1,bg_number] = y_pos
	    local_pos[2,bg_number] = z_pos
	    local_vel[0,bg_number] = 0.0
	    local_vel[1,bg_number] = 0.0
	    local_vel[2,bg_number] = 0.0
	    local_mass[ bg_number] = mean(m0)					; give them mean mass
	    local_rho[  bg_number] = mean(m0)*nbackground_base^3/(20000.0^3)	; give them reasonable density
	    local_u[    bg_number] = max(internal_energy_gas)			; give them maximum temperature
	    bg_number += 1
	  endif
        endfor
      endfor
    endfor

  local_pos = local_pos[*,0:bg_number-1]
  local_vel = local_vel[*,0:bg_number-1]
  local_mass= local_mass[0:bg_number-1]
  local_rho = local_rho[0:bg_number-1]
  local_u   = local_u[0:bg_number-1]

  
  for iii=0,2 do begin
    local_pos[iii,*] += mean(halo_position[iii])
    index = where(local_pos[iii,*] lt 0)
    if index[0] ne -1 then local_pos[iii,index] += boxsize
    index = where(local_pos[iii,*] gt boxsize)
    if index[0] ne -1 then local_pos[iii,index] -= boxsize
  endfor

endif

if n_elements(x0) gt 1 then N0 = n_elements(x0) + bg_number else n0 = 0
if n_elements(x1) gt 1 then N1 = n_elements(x1) else n1 = 0
if n_elements(x2) gt 1 then N2 = n_elements(x2) else n2 = 0
if n_elements(x3) gt 1 then N3 = n_elements(x3) else n3 = 0
if n_elements(x4) gt 1 then N4 = n_elements(x4) else n4 = 0
if n_elements(x5) gt 1 then N5 = n_elements(x5) else n5 = 0

print,' '
print,' '
print,n0,n1,n2,n3,n4,n5
print,' '
print,' '

if N0 eq 1 then N0= 0
if N1 eq 1 then N1= 0
if N2 eq 1 then N2= 0
if N3 eq 1 then N3= 0
if N4 eq 1 then N4= 0
if N_ELEMENTS(x5) gt 0 then begin
    if x5[0] eq 0 then N5= 0
endif



npart 	= [N0,N1,N2,N3,N4,N5]
npartall= [N0,N1,N2,N3,N4,N5]

N_total = N0 + N1 + N2 + N3 + N4 + N5


id 		= l64indgen(N_total)+1
pos 		= fltarr(3, N_total)
vel		= fltarr(3, N_total)
u 		= fltarr(N0)*0.0+1.0
rho 		= fltarr(N0)
mass 		= fltarr(N_total)

if keyword_set(rho_gas) 		then rho[0:n_elements(rho_gas)-1] = rho_gas
if keyword_set(internal_energy_gas) 	then u[  0:n_elements(rho_gas)-1] = internal_energy_gas

if N0 gt 0 then begin
    pos[0,0:n_elements(x0)-1] = x0
    pos[1,0:n_elements(x0)-1] = y0
    pos[2,0:n_elements(x0)-1] = z0
    vel[0,0:n_elements(x0)-1] = vx0
    vel[1,0:n_elements(x0)-1] = vy0
    vel[2,0:n_elements(x0)-1] = vz0
    mass[0:n_elements(x0)-1] = m0

    if not keyword_set(sunrise_only) then begin

      pos[0,n_elements(x0):N0-1] = local_pos[0,*]
      pos[1,n_elements(x0):N0-1] = local_pos[1,*]
      pos[2,n_elements(x0):N0-1] = local_pos[2,*]

      vel[0,n_elements(x0):N0-1] = local_vel[0,*]
      vel[1,n_elements(x0):N0-1] = local_vel[1,*]
      vel[2,n_elements(x0):N0-1] = local_vel[2,*]

      rho[ n_elements(x0):N0-1] = local_rho
      mass[n_elements(x0):N0-1] = local_mass
      u[   n_elements(x0):N0-1] = local_u
    endif
;local_mass= local_mass[0:bg_number-1]
;  local_rho = local_rho[0:bg_number-1]
;  local_u   = local_u[0:bg_number-1]

endif

if N1 gt 0 then begin
  pos[0,N0:N0+N1-1] = x1
  pos[1,N0:N0+N1-1] = y1
  pos[2,N0:N0+N1-1] = z1
  vel[0,N0:N0+N1-1] = vx1
  vel[1,N0:N0+N1-1] = vy1
  vel[2,N0:N0+N1-1] = vz1

  mass[N0:N0+N1-1] = m1
endif

if N2 gt 0 then begin
  pos[0,N0+N1:N0+N1+N2-1] = x2
  pos[1,N0+N1:N0+N1+N2-1] = y2
  pos[2,N0+N1:N0+N1+N2-1] = z2
  vel[0,N0+N1:N0+N1+N2-1] = vx2
  vel[1,N0+N1:N0+N1+N2-1] = vy2
  vel[2,N0+N1:N0+N1+N2-1] = vz2

  mass[N0+N1:N0+N1+N2-1] = m2
endif

if N3 gt 0 then begin
  pos[0,N0+N1+N2:N0+N1+N2+N3-1] = x3
  pos[1,N0+N1+N2:N0+N1+N2+N3-1] = y3
  pos[2,N0+N1+N2:N0+N1+N2+N3-1] = z3
  vel[0,N0+N1+N2:N0+N1+N2+N3-1] = vx3
  vel[1,N0+N1+N2:N0+N1+N2+N3-1] = vy3
  vel[2,N0+N1+N2:N0+N1+N2+N3-1] = vz3

  mass[N0+N1+N2:N0+N1+N2+N3-1] = m3
endif

if N4 gt 0 then begin
  pos[0,N0+N1+N2+N3:N0+N1+N2+N3+N4-1] = x4
  pos[1,N0+N1+N2+N3:N0+N1+N2+N3+N4-1] = y4
  pos[2,N0+N1+N2+N3:N0+N1+N2+N3+N4-1] = z4
  vel[0,N0+N1+N2+N3:N0+N1+N2+N3+N4-1] = vx4
  vel[1,N0+N1+N2+N3:N0+N1+N2+N3+N4-1] = vy4
  vel[2,N0+N1+N2+N3:N0+N1+N2+N3+N4-1] = vz4

  mass[N0+N1+N2+N3:N0+N1+N2+N3+N4-1] = m4
endif

if N5 gt 0 then begin
  pos[0,N0+N1+N2+N3+N4:N0+N1+N2+N3+N4+N5-1] = x5
  pos[1,N0+N1+N2+N3+N4:N0+N1+N2+N3+N4+N5-1] = y5
  pos[2,N0+N1+N2+N3+N4:N0+N1+N2+N3+N4+N5-1] = z5
  vel[0,N0+N1+N2+N3+N4:N0+N1+N2+N3+N4+N5-1] = vx5
  vel[1,N0+N1+N2+N3+N4:N0+N1+N2+N3+N4+N5-1] = vy5
  vel[2,N0+N1+N2+N3+N4:N0+N1+N2+N3+N4+N5-1] = vz5

  mass[N0+N1+N2+N3+N4:N0+N1+N2+N3+N4+N5-1] = m5
endif


true_snapshot_format = 1
;true_snapshot_format = 0


if keyword_set(true_snapshot_format) then begin
  filename = snapshot_name
  fid = h5f_create(filename)
  head = h5g_create(fid,'Header')

  datatype = h5t_idl_create(npart)
  npart_thisfile = h5a_create(head,'NumPart_ThisFile',datatype,h5s_create_simple(n_elements(npart)))
  h5a_write,npart_thisfile,npart
  h5a_close,npart_thisfile

  print,npart

  datatype = h5t_idl_create(npart)
  npart_total = h5a_create(head,'NumPart_Total',datatype,h5s_create_simple(n_elements(npart)))
  h5a_write,npart_total,npart
  h5a_close,npart_total

  ;print,npart

  datatype = h5t_idl_create(npart)
  npart_total = h5a_create(head,'NumPart_Total_HighWord',datatype,h5s_create_simple(n_elements(npart)))
  h5a_write,npart_total,npart*0
  h5a_close,npart_total

  ;print,npart*0

  datatype = h5t_idl_create(massarr)
  MassTable = h5a_create(head,'MassTable',datatype,h5s_create_simple(n_elements(massarr)))
  h5a_write,MassTable,massarr
  h5a_close,MassTable

  print,massarr

  datatype = h5t_idl_create(time)
  Time_ = h5a_create(head,'Time',datatype,h5s_create_simple(n_elements(time)))
  h5a_write,Time_,time
  h5a_close,Time_

  ;print,time

  datatype = h5t_idl_create(redshift)
  redshift_ = h5a_create(head,'Redshift',datatype,h5s_create_simple(n_elements(redshift)))
  h5a_write,redshift_,redshift
  h5a_close,redshift_

  ;print,redshift

  datatype = h5t_idl_create(boxsize)
  boxsize_ = h5a_create(head,'BoxSize',datatype,h5s_create_simple(n_elements(boxsize)))
  h5a_write,boxsize_,boxsize
  h5a_close,boxsize_

  ;print,boxsize

  numfilespersnapshot = 1
  datatype = h5t_idl_create(numfilespersnapshot)
  numfilespersnapshot_ = h5a_create(head,'NumFilesPerSnapshot',datatype,h5s_create_simple(n_elements(numfilespersnapshot)))
  h5a_write,numfilespersnapshot_,numfilespersnapshot
  h5a_close,numfilespersnapshot_

  ;print,numfilespersnapshot
  
  Omega0 = 0.27
  datatype = h5t_idl_create(Omega0)
  Omega0_ = h5a_create(head,'Omega0',datatype,h5s_create_simple(n_elements(Omega0)))
  h5a_write,Omega0_,Omega0
  h5a_close,Omega0_

  ;print,Omega0

  OmegaLambda = 0.73
  datatype = h5t_idl_create(OmegaLambda)
  OmegaLambda_ = h5a_create(head,'OmegaLambda',datatype,h5s_create_simple(n_elements(OmegaLambda)))
  h5a_write,OmegaLambda_,OmegaLambda
  h5a_close,OmegaLambda_

  ;print,OmegaLambda

  HubbleParam = 0.704000d
  datatype = h5t_idl_create(HubbleParam)
  HubbleParam_ = h5a_create(head,'HubbleParam',datatype,h5s_create_simple(n_elements(HubbleParam)))
  h5a_write,HubbleParam_,HubbleParam
  h5a_close,HubbleParam_

  ;print,HubbleParam

  Flag_SFR = 1
  datatype = h5t_idl_create(Flag_SFR)
  Flag_SFR_ = h5a_create(head,'Flag_Sfr',datatype,h5s_create_simple(n_elements(Flag_SFR)))
  h5a_write,Flag_SFR_,Flag_SFR
  h5a_close,Flag_SFR_

  Flag_Metals = 0
  datatype = h5t_idl_create(Flag_Metals)
  Flag_SFR_ = h5a_create(head,'Flag_Metals',datatype,h5s_create_simple(n_elements(Flag_SFR)))
  h5a_write,Flag_SFR_,Flag_SFR
  h5a_close,Flag_SFR_

  Flab_Cooling = 1
  datatype = h5t_idl_create(Flag_Cooling)
  Flag_Cooling_ = h5a_create(head,'Flag_Cooling',datatype,h5s_create_simple(n_elements(Flag_Cooling)))
  h5a_write,Flag_Cooling_,Flag_Cooling
  h5a_close,Flag_Cooling_

  Flag_StellarAge = 1
  datatype = h5t_idl_create(Flag_StellarAge)
  Flag_StellarAge_ = h5a_create(head,'Flag_StellarAge',datatype,h5s_create_simple(n_elements(Flag_StellarAge)))
  h5a_write,Flag_StellarAge_,Flag_StellarAge
  h5a_close,Flag_StellarAge_

  Flag_Feedback = 1
  datatype = h5t_idl_create(Flag_Feedback)
  Flag_Feedback_ = h5a_create(head,'Flag_Feedback',datatype,h5s_create_simple(n_elements(Flag_Feedback)))
  h5a_write,Flag_Feedback_,Flag_Feedback
  h5a_close,Flag_Feedback_

  Flag_DoublePrecision = 0
  datatype = h5t_idl_create(Flag_DoublePrecision)
  Flag_DoublePrecision_ = h5a_create(head,'Flag_DoublePrecision',datatype,h5s_create_simple(n_elements(Flag_DoublePrecision)))
  h5a_write,Flag_DoublePrecision_,Flag_DoublePrecision
  h5a_close,Flag_DoublePrecision_

  h5g_close,head


  for iii=0,5 do begin
      print, npart[iii]
    if npart[iii] gt 0 then begin
      groupname='PartType'+strcompress(string(iii),/remove_all)
      groupid = h5g_create(fid,groupname)

      if iii eq 0 then coordinates = transpose([[x0],[y0],[z0]])
      if iii eq 1 then coordinates = transpose([[x1],[y1],[z1]])
      if iii eq 2 then coordinates = transpose([[x2],[y2],[z2]])
      if iii eq 3 then coordinates = transpose([[x3],[y3],[z3]])
      if iii eq 4 then coordinates = transpose([[x4],[y4],[z4]])
      if iii eq 5 then coordinates = transpose([[x5],[y5],[z5]])

;print,' '
;print,' '
;print,' '
;help,x0
;help,x1
;help,x4
;print, coordinates[1,0:5]
;	help,coordinates
;print,' '
;print,' '
;print,' '


      ;write_x = coordinates[0,*]
      ;write_y = coordinates[1,*]
      ;write_z = coordinates[2,*]
      ;DISTANCE = SQRT(write_x^2 + write_y^2 + write_z^2)
      ;NUMBER = N_ELEMENTS(DISTANCE)
      ;UNIQUE_DISTANCE = UNIQ(DISTANCE,SORT(DISTANCE))
      ;NEWNUMBER = N_ELEMENTS(UNIQUE_DISTANCE)
      ;print, 'NUMBER - NEWNUMBER =   ', NUMBER-NEWNUMBER, '    PARTICLE TYPE  ', iii,   '    TEST   NUMBER=   ', NUMBER

      datatype = h5t_idl_create(coordinates)
      Coordinates_ = h5d_create(groupid,'Coordinates',datatype,h5s_create_simple([3,n_elements(Coordinates)/3]))
      h5d_write,Coordinates_,Coordinates
      h5d_close,Coordinates_

      if iii eq 0 then velocities = transpose([[vx0],[vy0],[vz0]])
      if iii eq 1 then velocities = transpose([[vx1],[vy1],[vz1]])
      if iii eq 2 then velocities = transpose([[vx2],[vy2],[vz2]])
      if iii eq 3 then velocities = transpose([[vx3],[vy3],[vz3]])
      if iii eq 4 then velocities = transpose([[vx4],[vy4],[vz4]])
      if iii eq 5 then velocities = transpose([[vx5],[vy5],[vz5]])

      datatype = h5t_idl_create(velocities)
      velocities_ = h5d_create(groupid,'Velocities',datatype,h5s_create_simple([3,n_elements(velocities)/3]))
      h5d_write,velocities_,velocities
      h5d_close,velocities_

      long_ids = 0
      if keyword_set(long_ids) then begin
        if iii eq 0 then ids = l64indgen(n_elements(x0))
        if iii eq 1 then ids = l64indgen(n_elements(x1)) + n_elements(x0) + 1
        if iii eq 2 then ids = l64indgen(n_elements(x2)) + n_elements(x0) + n_elements(x1) + 1
        if iii eq 3 then ids = l64indgen(n_elements(x3)) + n_elements(x0) + n_elements(x1) + n_elements(x2) +1
        if iii eq 4 then ids = l64indgen(n_elements(x4)) + n_elements(x0) + n_elements(x1) + n_elements(x2) + n_elements(x3) +1
        if iii eq 5 then ids = l64indgen(n_elements(x5)) + n_elements(x0) + n_elements(x1) + n_elements(x2) + n_elements(x3) + n_elements(x4) + 1
      endif else begin
        if iii eq 0 then ids = id0
        if iii eq 1 then ids = id1
        if iii eq 2 then ids = id2
        if iii eq 3 then ids = id3
        if iii eq 4 then ids = id4
        if iii eq 5 then ids = id5
      endelse

;help,ids

      datatype = h5t_idl_create(ids)
      ids_ = h5d_create(groupid,'ParticleIDs',datatype,h5s_create_simple(n_elements(ids)))
      h5d_write,ids_,ids
      h5d_close,ids_

      if massarr[iii] eq 0 then begin
        if iii eq 0 then masses = m0
        if iii eq 1 then masses = m1
        if iii eq 2 then masses = m2
        if iii eq 3 then masses = m3
        if iii eq 4 then masses = m4
        if iii eq 5 then masses = m5

        datatype = h5t_idl_create(masses)
        masses_ = h5d_create(groupid,'Masses',datatype,h5s_create_simple(n_elements(masses)))
        h5d_write,masses_,masses
        h5d_close,masses_
      endif

      if iii eq 0 and keyword_set(met0) then begin
        datatype = h5t_idl_create(met0)
        metallicity_ = h5d_create(groupid,'GFM_Metallicity',datatype,h5s_create_simple(n_elements(met0)))
        h5d_write,metallicity_,met0
        h5d_close,metallicity_
      endif

      if iii eq 4 and keyword_set(met4) then begin
        datatype = h5t_idl_create(met4)
        metallicity_ = h5d_create(groupid,'GFM_Metallicity',datatype,h5s_create_simple(n_elements(met4)))
        h5d_write,metallicity_,met4
        h5d_close,metallicity_
      endif




      if iii eq 4 and keyword_set(gfm_initial_mass) then begin
        datatype = h5t_idl_create(gfm_initial_mass)
        initial_mass_ = h5d_create(groupid,'GFM_InitialMass',datatype,h5s_create_simple(n_elements(gfm_initial_mass)))
        h5d_write,initial_mass_,gfm_initial_mass
        h5d_close,initial_mass_
      endif

    
      if iii eq 4 then begin
	if npart[4] gt 0 then begin
  	  datatype = h5t_idl_create(stellarage)
          stellarage_ = h5d_create(groupid,'GFM_StellarFormationTime',datatype,h5s_create_simple(n_elements(stellarage)))
          h5d_write,stellarage_,stellarage
          h5d_close,stellarage_
	endif
      endif


      if iii eq 0 then begin
	if npart[0] gt 0 then begin
if keyword_set(rho_gas) then begin
  	  datatype = h5t_idl_create(rho_gas)
          density_ = h5d_create(groupid,'Density',datatype,h5s_create_simple(n_elements(rho_gas)))
          h5d_write,density_,rho_gas
          h5d_close,density_
endif


  	  datatype = h5t_idl_create(internal_energy_gas)
          internalenergy_ = h5d_create(groupid,'InternalEnergy',datatype,h5s_create_simple(n_elements(internal_energy_gas)))
          h5d_write,internalenergy_,internal_energy_gas
          h5d_close,internalenergy_

if keyword_set(smoothinglength) then begin
  	  datatype = h5t_idl_create(smoothinglength)
          smoothinglength_ = h5d_create(groupid,'SmoothingLength',datatype,h5s_create_simple(n_elements(smoothinglength)))
          h5d_write,smoothinglength_,smoothinglength
          h5d_close,smoothinglength_
endif 

if keyword_set(sfr) then begin
  	  datatype = h5t_idl_create(sfr)
          smoothinglength_ = h5d_create(groupid,'StarFormationRate',datatype,h5s_create_simple(n_elements(sfr)))
          h5d_write,smoothinglength_,sfr
          h5d_close,smoothinglength_
endif

if keyword_set(ea) then begin
  	  datatype = h5t_idl_create(ea)
          smoothinglength_ = h5d_create(groupid,'ElectronAbundance',datatype,h5s_create_simple(n_elements(ea)))
          h5d_write,smoothinglength_,ea
          h5d_close,smoothinglength_
endif 

if keyword_set(volume) then begin
  	  datatype = h5t_idl_create(volume)
          smoothinglength_ = h5d_create(groupid,'Volume',datatype,h5s_create_simple(n_elements(volume)))
          h5d_write,smoothinglength_,volume
          h5d_close,smoothinglength_
endif

	endif
      endif



    endif
  endfor

  h5f_close,fid


endif else begin
  openw,1,fout,/f77_unformatted
  writeu,1, npart,massarr,time,redshift,flag_sfr,flag_feedback,npartall,flag_cooling,num_files,BoxSize,la
  writeu,1, pos
  writeu,1, vel

;  id = fix(id)
  writeu,1, id

  help,id
  writeu,1, mass
  writeu,1, u
  writeu,1, rho
  close,1
endelse





end


