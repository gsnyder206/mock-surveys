;####
;#### name:   fload_cylindrical_volume.pro
;#### author: Greg Snyder (gsnyder@stsci.edu)
;#### purpose: subroutine to load a cylindrical volume from a snapshot
;#### calls:  IDL HDF5 sub-routines
;####


function fload_cylindrical_volume, frun, num, parttype, BoxSize, RadSize, v_Ingress, v_Egress, mass=mass,time=time,redshift=redshift,$
			x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,rho=rho,u=u,ids=ids,$
			metallicity=metallicity,$
			stellarage=stellarage,$
			hsml=hsml


  g_minus_1= 5.0/3.0-1.0
  PROTONMASS = 1.6726e-24
  BoltzMann_ergs= 1.3806d-16
;  MeanWeight= 0.6*PROTONMASS
  xh = 0.76
  yh = (1.0-xh)/(4.0*xh)



    if not keyword_set(num) then num=0
    if not keyword_set(frun) then frun='./'

    exts='0000'+strcompress(string(num),/remove_all)
    if num ge 1000 then do_four= 1
    if not keyword_set(do_four) then begin
      exts=strmid(exts,strlen(exts)-3,3)
    endif else begin
      exts=strmid(exts,strlen(exts)-4,4)
    endelse

    basename=strcompress(frun+'snapdir_'+exts+'/snap_*',/remove_all)
    filelist = findfile(basename)

    if filelist[0] eq '' then begin
	print," "
	print," "
	print,"Can't find snapshots. "
	print,"Searching in: "+basename
	print," "
	print," "
	return,-1
    endif

    n_files = n_elements(filelist)
    file_number = intarr(n_files)
    for i=0,n_files-1 do begin
      file_number[i] =  round(double( strmid(filelist[i],strpos(filelist[i],'snap_')+9,strpos(filelist[i],'.hdf5')-strpos(filelist[i],'snap_')-9) ) )
    endfor

    filelist = filelist[sort(file_number)]

    fname=strcompress(filelist[0],/remove_all)
    file_id= h5f_open(fname)
    hdr_group_id= h5g_open(file_id,"Header")

    npartTotal	= h5a_read(h5a_open_name(hdr_group_id,"NumPart_Total"))
    time     	= h5a_read(h5a_open_name(hdr_group_id,"Time"))
    redshift   	= h5a_read(h5a_open_name(hdr_group_id,"Redshift"))
    massarr	= h5a_read(h5a_open_name(hdr_group_id,"MassTable"))
    h5g_close, hdr_group_id
    h5f_close, file_id

    if npartTotal[parttype] eq 0 then return,0

    count = 0L

;    id 	        = lonarr(npartTotal[parttype])
;    return_val1 = dblarr(npartTotal[parttype])
;    return_val2 = dblarr(npartTotal[parttype])
;    return_val3 = dblarr(npartTotal[parttype])
    mass 	= dblarr(npartTotal[parttype])

    id_all = 0
    x_all = 0
    y_all = 0
    z_all = 0
    vx_all = 0
    vy_all = 0
    vz_all = 0
    mass = 0

    ;print,"Loading cylindrical volume..."
    ;print, v_Ingress, v_Egress
    for ii=0,n_elements(filelist)-1 do begin
      fname=strcompress(filelist[ii],/remove_all)
      print, "Reading HDF5 file: ", fname

      file_id		= h5f_open(fname)
      hdr_group_id	= h5g_open(file_id,"Header")
      npart		= h5a_read(h5a_open_name(hdr_group_id,"NumPart_ThisFile"))
      numfles		= h5a_read(h5a_open_name(hdr_group_id,"NumFilesPerSnapshot"))
      fileboxsize       = h5a_read(h5a_open_name(hdr_group_id,"BoxSize"))
      filehubble        = h5a_read(h5a_open_name(hdr_group_id,"HubbleParam"))
      h5g_close, hdr_group_id

      if npart[parttype] eq 0 then continue

      GFS_ASSERT, filehubble, 0.704d, 'Hubble Weirdness'
      groupName 	= strcompress('PartType'+string(parttype),/remove_all)
      group_id 		= h5g_open(file_id,groupName)
      xyz		= h5d_read(h5d_open(group_id,"Coordinates"))

      ;xyz[0,*] -= mean(pos[0])  ;center coordinates on input position
      ;xyz[1,*] -= mean(pos[1])
      ;xyz[2,*] -= mean(pos[2])

                                ;make a full 20mpc^3 box centered on
                                ;target position (I think?)
                                ;"periodicalize"


                                ;a way to do this more generally is to
                                ;assume some size scale R for the region
                                ;we want to excise, then say if (0 < x
                                ;< R, x += L, etc) oh but this
                                ;*replaces* the particle, it doesn't
                                ;add a new one... will have to think
                                ;on that...

      ; what about... center on mid-point?? or one of the *gress points?
      ; maybe center on (x_mid,y_mid,0), and only replicate in x,y?
      ; thought that's less general...

      ; this might work if we enforce only the extraction of 
      ; a cylinder of primary axis height equal to the sep.
      ; between the *gress points in that axis

                                ; better still:  extract exactly 1
                                ; cylinder from each snap along the
                                ; viewing direction, of length =
                                ; distance through box along LOS
      ; in the latter case, will want to center on midpoint
      ; and then replicate in all directions, then cut at end

      x_new = REFORM(xyz[0,*])
      y_new = REFORM(xyz[1,*])
      z_new = REFORM(xyz[2,*])

      N_orig = N_ELEMENTS(xyz[0,*])

      size=RadSize
      L = BoxSize

      orig_indices = where(xyz[0,*] gt -1d12)
      ;print, N_ELEMENTS(orig_indices), N_orig, N_ELEMENTS(orig_indices[UNIQ(orig_indices,SORT(orig_indices))])
      i_new = REFORM(orig_indices)

      index = where(xyz[0,*] lt size or xyz[0,*] gt (L-size), count )
      if count gt 0 then begin
          sign = (-1.0d)*(xyz[0,index] - size)/ABS(xyz[0,index] - size)
          x_new = [x_new,REFORM(xyz[0,index]+sign*L)]
          y_new = [y_new,REFORM(xyz[1,index])]
          z_new = [z_new,REFORM(xyz[2,index])]
          i_new = [i_new,REFORM(orig_indices[index])] ; list of indices into xyz[0,*]
          ;print, count, MAX(REFORM(xyz[0,*])), MIN(REFORM(xyz[0,*])), MAX(x_new), MIN(x_new)
      endif

      index = where(xyz[1,*] lt size or xyz[1,*] gt (L-size), count )
      if count gt 0 then begin
          sign = (-1.0d)*(xyz[1,index] - size)/ABS(xyz[1,index] - size)
          x_new = [x_new,REFORM(xyz[0,index])]
          y_new = [y_new,REFORM(xyz[1,index]+sign*L)]
          z_new = [z_new,REFORM(xyz[2,index])]
          i_new = [i_new,REFORM(orig_indices[index])]
          ;print, count, MAX(REFORM(xyz[1,*])), MIN(REFORM(xyz[1,*])), MAX(y_new), MIN(y_new)
      endif

      index = where(xyz[2,*] lt size or xyz[2,*] gt (L-size), count )
      if count gt 0 then begin
          sign = (-1.0d)*(xyz[2,index] - size)/ABS(xyz[2,index] - size)
          x_new = [x_new,REFORM(xyz[0,index])]
          y_new = [y_new,REFORM(xyz[1,index])]
          z_new = [z_new,REFORM(xyz[2,index]+sign*L)]
          i_new = [i_new,REFORM(orig_indices[index])]
          ;print, count, MAX(REFORM(xyz[2,*])), MIN(REFORM(xyz[2,*])), MAX(z_new), MIN(z_new)
      endif

      ;print, N_ELEMENTS(i_new), N_ELEMENTS(i_new[UNIQ(i_new,SORT(i_new))]) ;~3/8 replicated, as expected

      GFS_ASSERT, L, fileboxsize
      ;help, i_new
      ;print, N_ELEMENTS(x_new), N_ELEMENTS(i_new), N_ELEMENTS(REFORM(xyz[0,*])), N_ELEMENTS(orig_indices), MAX(z_new), MIN(z_new), MAX(y_new), MIN(y_new), MAX(x_new), MIN(x_new)

      ; need way to ensure no ID replication... must test

      ;index = where(xyz[0,*] lt -10000.0)
      ;if index[0] ne -1 then xyz[0,index] += 20000.0
      ;index = where(xyz[1,*] lt -10000.0)
      ;if index[0] ne -1 then xyz[1,index] += 20000.0
      ;index = where(xyz[2,*] lt -10000.0)
      ;if index[0] ne -1 then xyz[2,index] += 20000.0

      ;index = where(xyz[0,*] gt 10000.0)
      ;if index[0] ne -1 then xyz[0,index] -= 20000.0
      ;index = where(xyz[1,*] gt 10000.0)
      ;if index[0] ne -1 then xyz[1,index] -= 20000.0
      ;index = where(xyz[2,*] gt 10000.0)
      ;if index[0] ne -1 then xyz[2,index] -= 20000.0


      ; this is the part we'd want to change
      ; HOWEVER... def. don't want to do a for loop... write an array-
      ; based cross product function

      V_VECTOR = v_Egress - v_Ingress
      U_x = v_Ingress[0] - x_new
      U_y = v_Ingress[1] - y_new
      U_z = v_Ingress[2] - z_new

      V_x = V_VECTOR[0]
      V_y = V_VECTOR[1]
      V_z = V_VECTOR[2]

      ; let's use squared values to let up on the compute time
      NormSquared = (V_x^2.0 + V_y^2.0 + V_z^2.0)
      
      b = V_x*U_x + V_y*U_y + V_z*U_z
      C_x = V_y*U_z - V_z*U_y
      C_y = V_z*U_x - V_x*U_z
      C_z = V_x*U_y - V_y*U_x

      dSqNum = (C_x^2.0 + C_y^2.0 + C_z^2.0)
      dSqDen = NormSquared
      
      tNum = (-1.0d)*b
      tDen = NormSquared

      ;help, dNum, dDen, tNum, tDen
      t_QUANT = tNum/tDen
      dSq_QUANT = dSqNum/dSqDen

      ;print, MAX(t_QUANT), MIN(t_QUANT), MAX(d_QUANT), MIN(d_QUANT)

      nearcyl_indices = where((dSq_QUANT lt size^2.0) and (t_QUANT lt 1.0d) and (t_QUANT gt 0.0d),countcyl)
      ;print, N_ELEMENTS(nearcyl_indices), N_ELEMENTS(nearcyl_indices[UNIQ(nearcyl_indices,SORT(nearcyl_indices))])
      ;print, countcyl
      ;continue

      ;r = sqrt(  xyz[0,*]^2 + xyz[1,*]^2 + xyz[2,*]^2  )
      ;local_index = where(r lt size)

      ;help,nearcyl_indices,i_new
      if nearcyl_indices[0] ne -1 then begin
        local_index = i_new[nearcyl_indices] ; local_indices is used into the original arrays below
        ;print, N_ELEMENTS(local_index),N_ELEMENTS(local_index[UNIQ(local_index,SORT(local_index))]) ; local_index has unexpected repeats for gas only!
        ; this is what we need for all non-coordinate quantities
	vxyz_local	= h5d_read(h5d_open(group_id,"Velocities"))
	if keyword_set(massarr[parttype]) then begin
	  mass_local = replicate(massarr[parttype],npart[parttype])
	endif else begin
	  mass_local = h5d_read(h5d_open(group_id,"Masses"))
	endelse

    	idtemp    = h5d_read(h5d_open(group_id,"ParticleIDs"))
	id_local	= idtemp[local_index]

        if parttype eq 0 then begin
	  rho_local = h5d_read(h5d_open(group_id,"Density"))
	  u_local   = h5d_read(h5d_open(group_id,"InternalEnergy"))
    	  nume_array= h5d_read(h5d_open(group_id,"ElectronAbundance"))
	  mu = (1.0+4.0*yh)/(1.0+yh*nume_array)
	  met	    = h5d_read(h5d_open(group_id,"GFM_Metallicity"))
	  ;met	    = h5d_read(h5d_open(group_id,"Metallicity"))
	  hsml_	    = h5d_read(h5d_open(group_id,"SmoothingLength"))

;			stellarage=stellarage,$

;	  u_local    =  mu[local_index]*protonmass/BoltzMann_ergs * g_minus_1 * u_local[local_index] * 1e+10
;
	  u_local    	= u_local[local_index]
	  rho_local  	= rho_local[local_index]
	  met_local 	= met[local_index]
	  hsml_local	= hsml_[local_index]
	  
	endif

	if parttype eq 4 then begin
	  stellarage_local = h5d_read(h5d_open(group_id,"GFM_StellarFormationTime"))
	  ;stellarage_local = h5d_read(h5d_open(group_id,"StellarFormationTime"))

	  stellarage_local = stellarage_local[local_index]

	  met	    = h5d_read(h5d_open(group_id,"GFM_Metallicity"))
	  ;met	    = h5d_read(h5d_open(group_id,"Metallicity"))

	  met_local 	= met[local_index]
	endif

        ;xyz_local  = xyz[*,local_index]
	vxyz_local = vxyz_local[*,local_index]
        mass_local = mass_local[local_index]


                                ; qua?  if there exist some particles
                                ; already loaded, add them to the
                                ; vectors?
                                ;aha, I get it.  we are looping
                                ;through the snapshot subfiles, so if
                                ;we've already started, then just add
                                ;everything else to the same vector to
                                ;be returned to the user

        if keyword_set(x) then begin  
	  x =   [x      ,x_new[nearcyl_indices] ];[x	,transpose(xyz_local[0,*]) + mean(pos[0])]
	  y =   [y      ,y_new[nearcyl_indices] ];[y	,transpose(xyz_local[1,*]) + mean(pos[1])] 
	  z =   [z      ,z_new[nearcyl_indices] ];[z	,transpose(xyz_local[2,*]) + mean(pos[2])]
	  vx =  [vx	,transpose(vxyz_local[0,*])]
	  vy =  [vy	,transpose(vxyz_local[1,*])]
	  vz =  [vz	,transpose(vxyz_local[2,*])]
	  ids=  [ids	,id_local]
	  mass = [mass	,mass_local]
	  if parttype eq 0 then begin
	    u    = [u,u_local]
	    rho  = [rho,rho_local]
	    hsml = [hsml,hsml_local]
	    metallicity = [metallicity,met_local]
	  endif
	  if parttype eq 4 then begin
	    stellarage = [stellarage,stellarage_local]
	    metallicity = [metallicity,met_local]
	  endif 
        endif else begin
	  x = x_new[nearcyl_indices];transpose(xyz_local[0,*]) + mean(pos[0])
	  y = y_new[nearcyl_indices];transpose(xyz_local[1,*]) + mean(pos[1])
	  z = z_new[nearcyl_indices];transpose(xyz_local[2,*]) + mean(pos[2])
	  vx = transpose(vxyz_local[0,*])
	  vy = transpose(vxyz_local[1,*])
	  vz = transpose(vxyz_local[2,*])
	  mass = mass_local
	  ids = id_local
	  if parttype eq 0 then begin
	    u = u_local
	    rho = rho_local
	    hsml = hsml_local
	    metallicity = met_local
	  endif
	  if parttype eq 4 then begin
	    stellarage = stellarage_local
	    metallicity = met_local
	  endif
        endelse
      endif	

      h5g_close, group_id
      ;help, metallicity
endfor


return,1
;return,return_val


end


