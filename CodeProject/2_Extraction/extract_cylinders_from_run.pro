

;####
;#### name:   extract_cylinders_from_run.pro
;#### author: Greg Snyder (gsnyder@stsci.edu)
;#### purpose: converts cubic hydro simulation into individual cylinders that subtend a lightcone.
;####        Lightcone definition format as in 1_Geometry.
;#### calls:  process_single_cylinder.pro
;####


pro extract_cylinders_from_run, lightcone_infofile, out_dir, run_dir=run_dir, BoxSize=BoxSize, fakeids=fakeids, save_cylinders=save_cylinders, fudgethez=fudgethez
astrolib

help, !CPU, /str

;gadget_run_dir =  '/n/hernquistfs1/mvogelsberger/ComparisonProject/512_20Mpc/Gadget/output/' 
;arepo_run_dir  =  '/n/hernquistfs1/mvogelsberger/ComparisonProject/512_20Mpc/Arepo_ENERGY/output/' 

if not keyword_set(run_dir) and keyword_set(gadget_run_dir) then run_dir=gadget_run_dir
if not keyword_set(run_dir) and keyword_set(arepo_run_dir) then run_dir=arepo_run_dir

print, 'Lightcone cylinder generation from file: '+lightcone_infofile
print, '   From run in directory: '+run_dir
print, '   Output will be placed in directory: '+out_dir

READCOL, lightcone_infofile, id, snaplabel, snapz, v_In_x, v_In_y, v_In_z, v_Out_x, v_Out_y, v_Out_z,$
  q,q,q,q,q,q,q,q,q,q,q,CylRad,FORMAT='I,A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D',/SILENT

Ncyl = N_ELEMENTS(CylRad)

print, "Number of Slices: ", N_ELEMENTS(CylRad)

;return


h=0.704d
if not keyword_set(BoxSize) then BoxSize=25000.0d
;RadSize=1.500d * 1000.0d * h ; units brain fart
;snapnum=314
;v_Ingress = [0.0d,0.0d,0.0d] * 1000.0d * h
;v_Egress = [2.0408163d,1.9047619d,28.57142857d] * 1000.0d * h

cwd = FILE_EXPAND_PATH(".")
;print, cwd
ensure_dir, out_dir
cd, out_dir
spawn, 'cp '+lightcone_infofile+' .'



textout_file = "snapshot_cylinder_extraction_output.txt"
OPENW, TEXT_OUT_UNIT, textout_file, /GET_LUN, WIDTH=10000, /APPEND
;!TEXTUNIT=TEXT_OUT_UNIT

;FORPRINT, data_vector, TEXTOUT=5, FORMAT="D"


printf, TEXT_OUT_UNIT, "## Extracted lightcone cylinders from run: "+run_dir

;fudgelist = [68,73,78,83,88,93,98,103,108,113,118,123,128,133,138,143,148,153,158,163,168,173,178,183,188,193,198,203,208]
;fudgelist = [40,50,54,58]


for s=0,Ncyl-1 do begin  ;for all needed slices
;for s=40,40 do begin  ; testing
if keyword_set(fudgelist) then begin
    ind = where(fudgelist eq s, countfudge)
    if countfudge ne 1 then continue
end


FLUSH, TEXT_OUT_UNIT



result = process_single_cylinder(s, run_dir, snaplabel, CylRad, v_In_x, v_In_y, v_In_z, v_Out_x, v_Out_y, v_Out_z, id, BoxSize, $
                                 Ngas=Ngas, Nstars=Nstars,N_neg_stars=N_neg_stars,N_nonunique=N_nonunique,N_corrected_nonunique=N_corrected_nonunique,RadSize=RadSize,filebase=filebase,snapnum=snapnum, fakeids=fakeids, save_cylinders=save_cylinders, fudgethez=fudgethez)


if result eq 2 then continue



printf, TEXT_OUT_UNIT, filebase+'.hdf5', snapnum, Ngas, Nstars, N_neg_stars, N_nonunique, N_corrected_nonunique, RadSize, FORMAT='(A25,I8,I20,4I15,F15.3)'

endfor                   ; end snapshot for loop


FREE_LUN, TEXT_OUT_UNIT

cd, cwd

end
