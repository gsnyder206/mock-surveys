
;####
;#### name:   fload_cylindrical_volume_by_subhalo.pro
;#### author: Greg Snyder (gsnyder@stsci.edu)
;#### purpose: subroutine to load a cylindrical volume from a snapshot
;#### calls:  Paul Torrey's cosmo-idl bitbucket code for GFM HDF5 snapshots (https://bitbucket.org/ptorrey/cosmo-idl)
;####


pro cylinder_process_subhalo, frun, num, parttype, partexist, shxyz_orig, newposvec, len_type_var, offset_type_var,x,y,z,vx,vy,vz,metallicity,mass,rho,u,hsml,stellarage,ids,ea,sfr,initmass

          ;help, partid
          ;if partid[0] eq 0 then continue ;;  I think this means no particles of the given type

          ;help, partid
          ;print, partid
          ;help, partid
          partid = fload_subhalo_particle_id(frun,num,len_type_var,offset_type_var)
          partxyz = fload_subhalo_particle_xyz(frun,num,len_type_var,offset_type_var)
          partvxyz = fload_subhalo_particle_vxyz(frun,num,len_type_var,offset_type_var)
          partmet = fload_subhalo_particle_metallicity(frun,num,len_type_var,offset_type_var)
          partmass = fload_subhalo_particle_mass(frun,num,len_type_var,offset_type_var)


      write_x = REFORM(partxyz[0,*])
      write_y = REFORM(partxyz[1,*])
      write_z = REFORM(partxyz[2,*])
      ;DISTANCE = SQRT(write_x^2 + write_y^2 + write_z^2)
      ;NUMBER = N_ELEMENTS(DISTANCE)
      ;UNIQUE_DISTANCE = UNIQ(DISTANCE,SORT(DISTANCE))
      ;NEWNUMBER = N_ELEMENTS(UNIQUE_DISTANCE)
      ;DBOOL = LONARR(NUMBER)+1L
      ;DBOOL[UNIQUE_DISTANCE]=0L
      ;ISAME = where(DBOOL eq 1L)

      ;same = where( ABS(DISTANCE - DISTANCE[ISAME[0]]) lt 0.00001d )

      ;print, '   '
      ;print, 'NUMBER - NEWNUMBER =   ', NUMBER-NEWNUMBER, '    PARTICLE TYPE  ', parttype
      ;help, len_type_var, offset_type_var
      ;print, len_type_var, offset_type_var
      ;print, 'INDICES   ', same
      ;print, write_x[same]
      ;print, write_y[same]
      ;print, write_z[same]
      ;print, 'NEWPOSVEC    ', newposvec
      ;print, 'ORIG SHXYZ    ', shxyz_orig
      ;print, '   '



          if parttype eq 0 then begin
              partrho = fload_subhalo_particle_rho(frun,num,len_type_var,offset_type_var)
              parttemp = fload_subhalo_particle_u(frun,num,len_type_var,offset_type_var)
              parthsml = fload_subhalo_particle_hsml(frun,num,len_type_var,offset_type_var)
              partea = fload_subhalo_particle_electron_abundance(frun,num,len_type_var,offset_type_var)
              partsfr = fload_subhalo_particle_sfr(frun,num,len_type_var,offset_type_var)
              partnha = fload_subhalo_particle_neutral_hydrogen_abundance(frun,num,len_type_var,offset_type_var)
          endif

          if parttype eq 4 then begin
              partage = fload_subhalo_particle_stellar_age(frun,num,len_type_var,offset_type_var)
              partinitmass = fload_subhalo_particle_stellar_initial_mass(frun,num,len_type_var,offset_type_var)
          endif




      if parttype eq 0 then begin
          if partexist eq 1 then begin
              rho = [TEMPORARY(rho),REFORM(partrho)]
              u = [TEMPORARY(u), REFORM(parttemp)]
              hsml = [TEMPORARY(hsml),REFORM(parthsml)]
              ea = [TEMPORARY(ea),REFORM(partea)]
              sfr = [TEMPORARY(sfr),REFORM(partsfr)]
              mass = [TEMPORARY(mass),REFORM(partmass*partnha)]    ; MP OFF version.  "MP ON" -->  (1 - partnha) ?

          endif else begin
              rho = REFORM(partrho)
              u = REFORM(parttemp)
              hsml = REFORM(parthsml)
              ea = REFORM(partea)
              sfr = REFORM(partsfr)
              mass = REFORM(partmass*partnha)

          endelse
      endif


      if parttype eq 4 then begin
          if partexist eq 1 then begin
              stellarage = [TEMPORARY(stellarage),REFORM(partage)]
              initmass = [TEMPORARY(initmass),REFORM(partinitmass)]
              mass = [TEMPORARY(mass),REFORM(partmass)]
          endif else begin
              stellarage = REFORM(partage)
              initmass = REFORM(partinitmass)
              mass = REFORM(partmass)
          endelse
      endif


      if partexist eq 1 then begin
          ids = [TEMPORARY(ids),REFORM(partid)]
          x = [TEMPORARY(x),REFORM(partxyz[0,*])-shxyz_orig[0]+newposvec[0] ] ;reposition
          y = [TEMPORARY(y),REFORM(partxyz[1,*])-shxyz_orig[1]+newposvec[1] ]
          z = [TEMPORARY(z),REFORM(partxyz[2,*])-shxyz_orig[2]+newposvec[2] ]
          vx = [TEMPORARY(vx),REFORM(partvxyz[0,*])]
          vy = [TEMPORARY(vy),REFORM(partvxyz[1,*])]
          vz = [TEMPORARY(vz),REFORM(partvxyz[2,*])]
          metallicity = [TEMPORARY(metallicity),REFORM(partmet)]


      endif else begin
          ids = REFORM(partid)
          x = REFORM(partxyz[0,*])-shxyz_orig[0]+newposvec[0] ;reposition
          y = REFORM(partxyz[1,*])-shxyz_orig[1]+newposvec[1]
          z = REFORM(partxyz[2,*])-shxyz_orig[2]+newposvec[2]
          vx = REFORM(partvxyz[0,*])
          vy = REFORM(partvxyz[1,*])
          vz = REFORM(partvxyz[2,*])
          metallicity = REFORM(partmet)
          partexist = 1
      endelse



      write_x = x
      write_y = y
      write_z = z
      ;DISTANCE = SQRT(write_x^2 + write_y^2 + write_z^2)
      ;NUMBER = N_ELEMENTS(DISTANCE)
      ;UNIQUE_DISTANCE = UNIQ(DISTANCE,SORT(DISTANCE))
      ;NEWNUMBER = N_ELEMENTS(UNIQUE_DISTANCE)
      ;DBOOL = LONARR(NUMBER)+1L
      ;DBOOL[UNIQUE_DISTANCE]=0L
      ;ISAME = where(DBOOL eq 1L)

      ;same = where( ABS(DISTANCE - DISTANCE[ISAME[0]]) lt 0.00001 )

      ;print, '   '
      ;print, 'NUMBER - NEWNUMBER =   ', NUMBER-NEWNUMBER, '    PARTICLE TYPE  ', parttype
      ;print, 'INDICES   ', same
      ;print, write_x[same]
      ;print, write_y[same]
      ;print, write_z[same]
      ;print, 'NEWPOSVEC    ', newposvec
      ;print, 'ORIG SHXYZ    ', shxyz_orig
      ;print, '   '




return
end







function fload_cylindrical_volume_by_subhalo, frun, num, parttype, BoxSize, RadSize, v_Ingress, v_Egress, mass=mass,time=time,redshift=redshift,$
			x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,rho=rho,u=u,ids=ids,$
			metallicity=metallicity,$
			stellarage=stellarage, $
			hsml=hsml, ea=ea, sfr=sfr, save_cylinders=save_cylinders, initmass=initmass, cylstring=cylstring


  g_minus_1= 5.0/3.0-1.0
  PROTONMASS = 1.6726e-24
  BoltzMann_ergs= 1.3806d-16
;  MeanWeight= 0.6*PROTONMASS
  xh = 0.76d
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

    ok = fload_subhalo_all(frun, num)

    len_type_subs       = fload_subhalo_len_type(1)
    offset_type_subs    = fload_subhalo_offset_type(1)
    shxyz                 = DOUBLE( fload_subhalo_pos(1) )

    shmstar             = fload_subhalo_stellar_mass(1)
    shmhalo             = fload_subhalo_mass(1)
    shZgas              = fload_subhalo_gas_metallicity(1)
    shZstar             = fload_subhalo_stellar_metallicity(1)
    shsfr               = fload_subhalo_sfr(1)
    shb                 = fload_subhalo_luminosity(1,/b_band)
    sh_r                = fload_subhalo_luminosity(1,/r_band)

    shgrnr              = fload_subhalo_grnr(1)
    grmass              = fload_group_mass(1)
    grm200              = fload_group_m200(1)
    grpos               = fload_group_pos(1)
    ;shb_dust            = fload_subhalo_luminosity(1, /b_band, /manual_extract, /dust_correction)


    mgi = where( ALOG10(grm200*(1.0e10)/0.704) gt 13.7, countmassivegroups)
    if countmassivegroups gt 1 then begin
        mgi_x = grpos[0,mgi]
        mgi_y = grpos[1,mgi]
        mgi_z = grpos[2,mgi]
    endif


    above_e12 = where(  ALOG10(grm200*(1.0e10)/0.704) gt 12.0, counte12  )

    print, countmassivegroups

    sh_groupm200 = grm200[shgrnr]

    help, shxyz

    if filelist[0] eq '' then begin
	print," "
	print," "
	print,"Can't find snapshots. "
	print,"Searching in: "+basename
	print," "
	print," "
	return,-1
    endif

    ;n_files = n_elements(filelist)
    ;file_number = intarr(n_files)
    ;for i=0,n_files-1 do begin
    ;  file_number[i] =  round(double( strmid(filelist[i],strpos(filelist[i],'snap_')+9,strpos(filelist[i],'.hdf5')-strpos(filelist[i],'snap_')-9) ) )
    ;endfor

    ;filelist = filelist[sort(file_number)]

    fname=strcompress(filelist[0],/remove_all)
    file_id= h5f_open(fname)
    hdr_group_id= h5g_open(file_id,"Header")

    npartTotal	= h5a_read(h5a_open_name(hdr_group_id,"NumPart_Total"))
    time     	= h5a_read(h5a_open_name(hdr_group_id,"Time"))
    redshift   	= h5a_read(h5a_open_name(hdr_group_id,"Redshift"))
    massarr	= h5a_read(h5a_open_name(hdr_group_id,"MassTable"))
    fileboxsize       = h5a_read(h5a_open_name(hdr_group_id,"BoxSize"))
    filehubble        = h5a_read(h5a_open_name(hdr_group_id,"HubbleParam"))
    GFS_ASSERT, filehubble, 0.704d, 'Hubble Weirdness'


    h5g_close, hdr_group_id
    h5f_close, file_id

    if npartTotal[parttype] eq 0 then return,0

    count = 0L

;    id 	        = lonarr(npartTotal[parttype])
;    return_val1 = dblarr(npartTotal[parttype])
;    return_val2 = dblarr(npartTotal[parttype])
;    return_val3 = dblarr(npartTotal[parttype])
    ;mass 	= dblarr(npartTotal[parttype])

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
    ;for ii=0,n_elements(filelist)-1 do begin
      ;fname=strcompress(filelist[ii],/remove_all)
      ;print, "Reading HDF5 file: ", fname

      ;file_id		= h5f_open(fname)
      ;hdr_group_id	= h5g_open(file_id,"Header")
      ;npart		= h5a_read(h5a_open_name(hdr_group_id,"NumPart_ThisFile"))
      ;numfles		= h5a_read(h5a_open_name(hdr_group_id,"NumFilesPerSnapshot"))
      ;fileboxsize       = h5a_read(h5a_open_name(hdr_group_id,"BoxSize"))
      ;filehubble        = h5a_read(h5a_open_name(hdr_group_id,"HubbleParam"))
      ;h5g_close, hdr_group_id

      ;if npart[parttype] eq 0 then continue

      ;groupName 	= strcompress('PartType'+string(parttype),/remove_all)
      ;group_id 		= h5g_open(file_id,groupName)
      ;xyz		= h5d_read(h5d_open(group_id,"Coordinates"))

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

      x_new = REFORM(shxyz[0,*])
      y_new = REFORM(shxyz[1,*])
      z_new = REFORM(shxyz[2,*])

      N_orig = N_ELEMENTS(shxyz[0,*])

      size=RadSize*1.0d
      L = BoxSize*1.0d

      orig_indices = LONG64( where(shxyz[0,*] gt -1d12, /L64) )
      ;print, N_ELEMENTS(orig_indices), N_orig, N_ELEMENTS(orig_indices[UNIQ(orig_indices,SORT(orig_indices))])
      i_new = REFORM(orig_indices)

      index = where(shxyz[0,*] lt size or shxyz[0,*] gt (L-size), count )
      if count gt 0 then begin
          sign = (-1.0d)*(shxyz[0,index] - size)/ABS(shxyz[0,index] - size)
          x_new = [x_new,REFORM(shxyz[0,index]+sign*L)]
          y_new = [y_new,REFORM(shxyz[1,index])]
          z_new = [z_new,REFORM(shxyz[2,index])]
          i_new = [i_new,REFORM(orig_indices[index])] ; list of indices into xyz[0,*]
          ;print, count, MAX(REFORM(xyz[0,*])), MIN(REFORM(xyz[0,*])), MAX(x_new), MIN(x_new)
      endif

      index = where(shxyz[1,*] lt size or shxyz[1,*] gt (L-size), count )
      if count gt 0 then begin
          sign = (-1.0d)*(shxyz[1,index] - size)/ABS(shxyz[1,index] - size)
          x_new = [x_new,REFORM(shxyz[0,index])]
          y_new = [y_new,REFORM(shxyz[1,index]+sign*L)]
          z_new = [z_new,REFORM(shxyz[2,index])]
          i_new = [i_new,REFORM(orig_indices[index])]
          ;print, count, MAX(REFORM(xyz[1,*])), MIN(REFORM(xyz[1,*])), MAX(y_new), MIN(y_new)
      endif

      index = where(shxyz[2,*] lt size or shxyz[2,*] gt (L-size), count )
      if count gt 0 then begin
          sign = (-1.0d)*(shxyz[2,index] - size)/ABS(shxyz[2,index] - size)
          x_new = [x_new,REFORM(shxyz[0,index])]
          y_new = [y_new,REFORM(shxyz[1,index])]
          z_new = [z_new,REFORM(shxyz[2,index]+sign*L)]
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
                                ;here, we should have the indices of
                                ;all the subhalos falling near the cylinder

      ;actually, I think these are indices into i_new

      ;help, i_new
      ;help, nearcyl_indices
      print, 'Number of subhalos: ', countcyl, N_orig
      if countcyl eq 0 then begin
          clear_common_blocks
          return,1
      endif

      partmask = LON64ARR(6)
      partmask[parttype]=1L
      ;print, partmask
      ;help, len_type_subs
      ;help, offset_type_subs

      ;x = DBLARR(100000000) ;& UNDEFINE(x)
      ;y = DBLARR(100000000) ;& UNDEFINE(y)
      ;z = DBLARR(100000000) ;& UNDEFINE(z)

      ;vx = DBLARR(100000000) ;& UNDEFINE(vx)
      ;vy = DBLARR(100000000) ;& UNDEFINE(vy)
      ;vz = DBLARR(100000000) ;& UNDEFINE(vz)

      ;metallicity = DBLARR(100000000) ;& UNDEFINE(metallicity)
      ;mass = DBLARR(100000000) ;& UNDEFINE(mass)
      ;ids = LONARR(100000000) ;& UNDEFINE(ids)

      ;rho = DBLARR(100000000) ;& UNDEFINE(rho)
      ;u = DBLARR(100000000) ;& UNDEFINE(u)
      ;hsml = DBLARR(100000000) ;& UNDEFINE(hsml)

      ;stellarage = DBLARR(100000000) ;& UNDEFINE(stellarage)

      ;partid = LONARR(100000000)
      ;partxyz = DBLARR(3,100000000)
      ;partvxyz = DBLARR(3,100000000)
      ;partmass = DBLARR(100000000)
      ;partrho = DBLARR(100000000)
      ;parttemp = DBLARR(100000000)
      ;parthsml = DBLARR(100000000)
      ;partage = DBLARR(100000000)




      shcat_file = 'cylinder_'+cylstring+'_subhalo_catalog.txt'
      OPENW, SH_OUT_UNIT, shcat_file, /GET_LUN, WIDTH=10000
      printf, SH_OUT_UNIT, "###  ", frun
      printf, SH_OUT_UNIT, '###  snapnum, shindex, D_KPC, x_0, y_0, z_0, x_new, y_new, z_new, Mstar, Mhalo, Zgas, Zstar, SFR, Bband, t_num, t_den, groupM200, isneargroup, mingroupdist, massive10Mpc, massive3Mpc, rband '


      shcount = 0L


      partexist = 0

      for i=0L,countcyl-1 do begin

          shi = i_new[nearcyl_indices[i]]
          shxyz_orig = REFORM(shxyz[*,shi])
          ;print, shi, nearcyl_indices[i], i

          ;print, len_type_subs[*,shi]
          if len_type_subs[parttype,shi] eq 0L then continue


          ;print out some subhalo lightcone catalog thing

          newposvec = [x_new[nearcyl_indices[i]] , y_new[nearcyl_indices[i]], z_new[nearcyl_indices[i]]]

                                ;snapshot, subhalo_index, xyz_orig, xyz_new, RA, DEC (x,
                                ;y?), distance (or "z" ?), vx, vy,
                                ;vz?,
                                ;M*, Mh, SFR, Z, Photometry?,
                                ;
          D_MEASURE_KPC = SQRT( dSq_QUANT[nearcyl_indices[i]])
          t_num_shi = tNum[nearcyl_indices[i]]
          t_den_shi = tDen[0] ;tDen[nearcyl_indices[i]]

          isneargroup = 0L
          mingroupdist = -1.0d
          if countmassivegroups gt 1 then begin
              D_list = SQRT( (shxyz_orig[0] - mgi_x)^2 + (shxyz_orig[1] - mgi_y)^2 + (shxyz_orig[2] - mgi_z)^2 )
              mingroupdist = MIN(D_list)
              if mingroupdist lt 10000.0 then isneargroup = 1L
          endif


          D_allgroups = SQRT( (shxyz_orig[0] - grpos[0,above_e12])^2 + (shxyz_orig[1] - grpos[1,above_e12])^2 + (shxyz_orig[2] - grpos[2,above_e12])^2 )
          nearby_gi = where(D_allgroups lt 10000.0d)
          close_gi = where(D_allgroups lt 3000.0d)

          mostmassive = -1.0
          mostmassive_close = -1.0

          mostmassive = MAX(grm200[above_e12[nearby_gi]])
          mostmassive_close = MAX(grm200[above_e12[close_gi]])

          if shmstar[shi] gt 1.0e-4 then printf, SH_OUT_UNIT, num, shi, D_MEASURE_KPC, shxyz_orig[0], shxyz_orig[1], shxyz_orig[2], newposvec[0], newposvec[1], newposvec[2], shmstar[shi], shmhalo[shi], shZgas[shi], shZstar[shi], shsfr[shi], shb[shi], t_num_shi, t_den_shi, sh_groupm200[shi], isneargroup, mingroupdist, mostmassive, mostmassive_close, sh_r[shi]



          if keyword_set(save_cylinders) then begin
          cylinder_process_subhalo, frun, num, parttype, partexist, shxyz_orig, newposvec, len_type_subs[*,shi]*partmask, offset_type_subs[*,shi],x,y,z,vx,vy,vz,metallicity,mass,rho,u,hsml,stellarage, ids, ea, sfr, initmass

          endif                     ; only if saving particles!

      ;if shcount mod 100 eq 0 then begin
          print, shcount
          ;help, x
          ;help, x, /MEMORY
          
          ;print, shxyz_orig, N_ELEMENTS(x)
      write_x = x
      write_y = y
      write_z = z
      ;DISTANCE = SQRT(write_x^2 + write_y^2 + write_z^2)
      ;NUMBER = N_ELEMENTS(DISTANCE)
      ;UNIQUE_DISTANCE = UNIQ(DISTANCE,SORT(DISTANCE))
      ;NEWNUMBER = N_ELEMENTS(UNIQUE_DISTANCE)
      ;DBOOL = LONARR(NUMBER)+1L
      ;DBOOL[UNIQUE_DISTANCE]=0L
      ;ISAME = where(DBOOL eq 1L)

      ;same = where( ABS(DISTANCE - DISTANCE[ISAME[0]]) lt 0.00001 )

      ;print, '   '
      ;print, 'NUMBER - NEWNUMBER =   ', NUMBER-NEWNUMBER, '    PARTICLE TYPE  ', parttype
      ;print, 'INDICES   ', same
      ;print, write_x[same]
      ;print, write_y[same]
      ;print, write_z[same]
      ;print, 'NEWPOSVEC    ', newposvec
      ;print, 'ORIG SHXYZ    ', shxyz_orig
      ;print, '   '

      ;endif

      ;boundary_check, x, y, z, boxsize  ; ??

          ;help, partxyz


      shcount = shcount + 1L
  endfor

FREE_LUN, SH_OUT_UNIT



      ;UNDEFINE, partxyz
      ;UNDEFINE, partvxyz
      ;UNDEFINE, partid
      ;UNDEFINE, partmet
      ;UNDEFINE, partmass
      ;UNDEFINE, partrho
      ;UNDEFINE, parttemp
      ;UNDEFINE, parthsml
      ;UNDEFINE, partage

  clear_common_blocks

  ;help
  return,1

end


