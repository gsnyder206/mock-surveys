
;####
;#### name:   process_single_cylinder.pro
;#### author: Greg Snyder (gsnyder@stsci.edu)
;#### purpose: subroutine to extract 1 cylinder from 1 cubic snapshot
;#### calls:  fload_cylindrical_volume_by_subhalo.pro  (calls cosmo-idl libraries, fload routines, may not be publicly available)
;####         OR  fload_cylindrical_volume.pro         (calls idl HDF5 routines directly)
;####


function process_single_cylinder, s, run_dir, snaplabel, CylRad, v_In_x, v_In_y, v_In_z, v_Out_x, v_Out_y, v_Out_z, id, BoxSize, $
                                 Ngas=Ngas, Nstars=Nstars,N_neg_stars=N_neg_stars,N_nonunique=N_nonunique,N_corrected_nonunique=N_corrected_nonunique,RadSize=RadSize,filebase=filebase,snapnum=snapnum, fakeids=fakeids, save_cylinders=save_cylinders, fudgethez=fudgethez




xyz  = ptrarr(3,6,/allocate_heap)
vxyz = ptrarr(3,6,/allocate_heap)
mass = ptrarr(6,/allocate_heap)
ids = ptrarr(6,/allocate_heap)

snapnum = LONG(STRTRIM(STRMID(snaplabel[s],1,3),2))


RadSize = CylRad[s]   ;*(1.01/1.1) ; NOTE: KLUGE!!!
v_Ingress = [v_In_x[s],v_In_y[s],v_In_z[s]]
v_Egress = [v_Out_x[s],v_Out_y[s],v_Out_z[s]]
print, "Snapshot number: "+STRING(snapnum)+"  with Ingress: "+STRING(v_Ingress[0])+" "+STRING(v_Ingress[1])+" "+STRING(v_Ingress[2])+"       and Egress: "+STRING(v_Egress[0])+" "+STRING(v_Egress[1])+" "+STRING(v_Egress[2])
filebase = 'cylinder_'+STRING(id[s],FORMAT='(I04)') ; literally matches with IDs from input file!
string_id_cylinder = STRING(id[s],FORMAT='(I04)')
result = file_info(filebase+'.hdf5')
if (result.EXISTS and keyword_set(save_cylinders)) then begin
    print, filebase+'.hdf5 exists, skipping...'
    return,2
endif
;continue
;printf, TEXT_OUT_UNIT, filebase+'.hdf5'
;continue



Ngas=0
Nstars=0
N_nonunique = 0
N_corrected_nonunique = 0

N_nu_gas = 0
N_nu_stars = 0
N_neg_stars = 0
for kk=0,5 do begin		; for all particle species
    ;continue
    if kk ne 0 and kk ne 4 then continue   ; and kk ne 1
    print, 'Particle type: ', kk
	  x = 0d
	  y = 0d
	  z = 0d
	  vx= 0d
	  vy= 0d
	  vz= 0d
	  mass_=0d
	  rho = 0d
	  u = 0d
          ;met4 = 0d
          ;met0 = 0d
          ;ages = 0d
          countnon=0L
          ids_ = 0L
          stellarage = 0d
          hsml = 0d
          ea = 0d
          sfr=0d
          initmass=0d


          if (not keyword_set(save_cylinders)) and (kk ne 0) then continue

	  ;pos = [SubhaloPos[0,k],SubhaloPos[1,k],SubhaloPos[2,k]]
	  xyz_temp    = fload_cylindrical_volume_by_subhalo(run_dir,snapnum,kk,BoxSize,RadSize,v_Ingress,v_Egress, $
				x=x,y=y,z=z,vx=vx,vy=vy,vz=vz, $
				mass=mass_,time=time_,redshift=redshift_, $
				rho=rho,u=u,ids=ids_, $
                                metallicity=metallicity, stellarage=stellarage, hsml=hsml, ea=ea, sfr=sfr, save_cylinders=save_cylinders, initmass=initmass, cylstring=string_id_cylinder)



          if keyword_set(save_cylinders) then begin   ; only if saving cylinders
          ;if kk eq 0 then met0=metallicity
          ;if kk eq 4 then met4=metallicity
          ;if kk eq 4 then ages=stellarage
          ;if kk eq 0 then gas_hsml=hsml

          ;if kk eq 0 then help,met0

	  ;index = where(x gt 20000.0)
	  ;if index[0] ne -1 then x[index] -= 20000.0
	  ;index = where(y gt 20000.0)
	  ;if index[0] ne -1 then y[index] -= 20000.0
	  ;index = where(z gt 20000.0)
	  ;if index[0] ne -1 then z[index] -= 20000.0

	  ;index = where(x lt 0.0)
	  ;if index[0] ne -1 then x[index] += 20000.0
	  ;index = where(y lt 0.0)
	  ;if index[0] ne -1 then y[index] += 20000.0
	  ;index = where(z lt 0.0)
	  ;if index[0] ne -1 then z[index] += 20000.0

          if n_elements(x) gt 1L then begin

            writeind = LINDGEN(n_elements(x))
            if kk eq 4 then begin
              allages = stellarage
              ;help, x
              ;help, allages
              pos_forma_ind = where(allages gt 0.0,countpos,NCOMPLEMENT=countneg)
                                ;with GFM on, "wind" particles are
                                ;also Type 4.  These record their
                                ;formation time as a negative
                                ;countdown-type thing, so their ages are
                                ;negatively-valued.  One day we want
                                ;to keep this information somehow as a
                                ;gas-like particle.  For now, we will
                                ;discard them.
              ;help, pos_forma_ind
              ;print, countpos
              if countpos gt 0L then begin
                  met4 = metallicity[pos_forma_ind]
                  smallmets = where(met4 lt 1.0d-10,countsmallmets) ; impose a metallicity floor for Sunrise usage?
                  if countsmallmets gt 0L then met4[smallmets] = 1.0d-10
                  Nstars = N_ELEMENTS(allages[pos_forma_ind])
                  ages = allages[pos_forma_ind]
                  writeind = pos_forma_ind
                  initmass = initmass[writeind]
              endif else continue
              N_neg_stars = countneg
              ;N_nu_stars = N_nonunique
            endif


                                ;deal with gas and/or star particles
                                ;that are repeated in the cylinder

                                ; this is a bit weird, since now we're
                                ; repeating volume in the same slice

                                ; however, number of repeats seems
                                ; small-ish during the galaxy era
                                ; (also mostly gas)

            writeids = ids_[writeind]

            if keyword_set(fakeids) then writeids = writeids*0L + LINDGEN(N_ELEMENTS(writeids))

            sind = SORT(writeids, /L64)
            sorted_ids = writeids[sind]
            unique_ind = UNIQ(sorted_ids)
            unique_ids = sorted_ids[unique_ind]
            bool = LONARR(N_ELEMENTS(sorted_ids)) + 1L
            bool[unique_ind]=0L
            nonunique_ind = where(bool eq 1L, countnon) ;indexes into sind

            N_nonunique =  N_nonunique + N_ELEMENTS(writeids) - N_ELEMENTS(writeids[UNIQ(writeids,SORT(writeids))])

            if countnon gt 0 then writeids[sind[nonunique_ind]] = MAX(writeids) + LINDGEN(N_ELEMENTS(nonunique_ind)) + 1L
            ;print, countnon, N_ELEMENTS(writeids) - N_ELEMENTS(writeids[UNIQ(writeids,SORT(writeids))])

            N_corrected_nonunique = N_corrected_nonunique + N_ELEMENTS(writeids) - N_ELEMENTS(writeids[UNIQ(writeids,SORT(writeids))])

            ;help, x
            ;help, writeind
            ;print, max(writeind)

            zfudge = 0.0d
            ;if (keyword_set(fudgethez) and kk eq 4) then zfudge = (0.1d )*RANDOMU(seed,N_ELEMENTS(writeind))


            DISTANCE = SQRT(x[writeind]^2 + y[writeind]^2 + z[writeind]^2 )
            UNIQUE_POS_INDICES = UNIQ(DISTANCE,SORT(DISTANCE)) ; indices into WRITEIND
            UBOOL = LON64ARR(N_ELEMENTS(writeind)) + 1L
            UBOOL[UNIQUE_POS_INDICES] = 0L
            help, DISTANCE, UNIQUE_POS_INDICES
            FIXIND = where(UBOOL eq 1L, COUNTFIX)
            ;if COUNTFIX gt 0 then zfudge = 0.5d + DBLARR(COUNTFIX)
            if COUNTFIX gt 0 then z[writeind[FIXIND]] = z[writeind[FIXIND]] + ( (0.5d )*RANDOMU(seed, COUNTFIX) - 0.25d)

            ;print, writeind[-1], writeind[-2]
            ;print, x[writeind[-1]], x[writeind[-2]]
	    *xyz[0,kk]  = x[writeind]
	    *xyz[1,kk]  = y[writeind]
	    *xyz[2,kk]  = z[writeind]
	    *vxyz[0,kk] = vx[writeind]
	    *vxyz[1,kk] = vy[writeind]
	    *vxyz[2,kk] = vz[writeind]
	    *mass[kk]   = mass_[writeind]
	    *ids[kk]   = writeids ;ids_[writeind]
            

            help, *xyz[2,kk]

            ;help, writeids

            ;print, (*ids[kk])[0:50]

            ;if countnon gt 0 then begin
            ;    print, (x[sind])[where(sorted_ids eq sorted_ids[nonunique_ind[5]])], (y[sind])[where(sorted_ids eq sorted_ids[nonunique_ind[5]])], (z[sind])[where(sorted_ids eq sorted_ids[nonunique_ind[5]])]
            ;endif

            

            if kk eq 0 then begin
              t_gas       = u[writeind]
              rho_gas     = rho[writeind]
              m_gas       = mass_[writeind]
              potential   = (mass_*0)[writeind]
              gas_hsml    = hsml[writeind]
              ea_gas      = ea[writeind]
              sfr_gas     = sfr[writeind]
              met0        = metallicity[writeind]
              Ngas        = N_ELEMENTS(rho_gas)
              ;N_nu_gas    = N_nonunique
            endif



	  endif



      endif   ; check if saving cylinders


;	  xyz_temp    = fload_limited_xyz(master_run_dir,snapnum,kk,'x',*SubhaloIDs[k],mass=mass_,time=time_,redshift=redshift_)
;	  if n_elements(xyz_temp) gt 3 then begin
;	    *xyz[0,kk]  = transpose(xyz_temp[0,*])
;	    *xyz[1,kk]  = transpose(xyz_temp[1,*])
;	    *xyz[2,kk]  = transpose(xyz_temp[2,*])
;	    *mass[kk]   = mass_
;	  endif

;	  xyz_temp   = fload_limited_vxyz(master_run_dir,snapnum,kk,'x',*SubhaloIDs[k])
;	  if n_elements(xyz_temp) gt 3 then begin
;	    *vxyz[0,kk] = transpose(xyz_temp[0,*])
;	    *vxyz[1,kk] = transpose(xyz_temp[1,*])
;	    *vxyz[2,kk] = transpose(xyz_temp[2,*])
;	  endif

endfor		; for the 6 particle types


;help, x, ages, met0, met4
if keyword_set(save_cylinders) then begin
if Ngas gt 0 and Nstars gt 0 then begin
internal_energy_gas = t_gas
pos = [0.0d,0.0d,0.0d]
repackage_mini_snapshot, $
  *xyz[0,0],*xyz[1,0],*xyz[2,0],*vxyz[0,0],*vxyz[1,0],*vxyz[2,0],*mass[0],*ids[0],$
  *xyz[0,1],*xyz[1,1],*xyz[2,1],*vxyz[0,1],*vxyz[1,1],*vxyz[2,1],*mass[1],*ids[1],$
  *xyz[0,2],*xyz[1,2],*xyz[2,2],*vxyz[0,2],*vxyz[1,2],*vxyz[2,2],*mass[2],*ids[2],$
  *xyz[0,3],*xyz[1,3],*xyz[2,3],*vxyz[0,3],*vxyz[1,3],*vxyz[2,3],*mass[3],*ids[3],$
  *xyz[0,4],*xyz[1,4],*xyz[2,4],*vxyz[0,4],*vxyz[1,4],*vxyz[2,4],*mass[4],*ids[4],$
  *xyz[0,5],*xyz[1,5],*xyz[2,5],*vxyz[0,5],*vxyz[1,5],*vxyz[2,5],*mass[5],*ids[5],$
  time_=time_,redshift_=redshift_,$
  internal_energy_gas=internal_energy_gas,$
  rho_gas = rho_gas,$
  filename=filebase+'.dat',$
  halo_position = pos, /true_snapshot_format, met0=met0, met4=met4, stellarage=ages, smoothinglength=gas_hsml, $
  snapshot_name=filebase+'.hdf5', $
  /sunrise_only, boxsize=BoxSize, ea=ea_gas, sfr=sfr_gas, initmass=initmass

print, 'Wrote: ', filebase+'.hdf5'
endif
endif ; only if saving cylinder snapshots



;do I want to free xyz, vxyz, mass, ids  here???
PTR_FREE,xyz,vxyz,mass,ids ;???
UNDEFINE,xyz
UNDEFINE,vxyz
UNDEFINE,mass
UNDEFINE,ids
UNDEFINE,x & UNDEFINE,y & UNDEFINE,z & UNDEFINE,vx & UNDEFINE,vy & UNDEFINE,vz & UNDEFINE,mass_ & UNDEFINE,ids_ & UNDEFINE, rho_gas & UNDEFINE, met0
UNDEFINE, met4 & UNDEFINE, ages & UNDEFINE, gas_hsml & UNDEFINE, internal_energy_gas & UNDEFINE, rho & UNDEFINE, hsml & UNDEFINE, stellarage & UNDEFINE, u & UNDEFINE, writeids & UNDEFINE, ea & UNDEFINE, sfr

HEAP_GC, /VERBOSE









return, 1
end
