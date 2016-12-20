
;#### cylinderscript*.pro
;#### example usage of extract_cylinders_from_run.pro



extract_cylinders_from_run, '~/cosmo_stuff/NatureLetter/BIG_UDF_FILES/hudf_75Mpc_11_10_xyz_AllSnaps.txt', '/n/ghernquist/gsnyder/cosmo/lightcones/Production/Illustris/L75n1820FP/hudf_75Mpc_11_10_xyz_COLDGAS', run_dir='/n/ghernquist/Illustris/Runs/L75n1820FP/', BoxSize=75000.0d, /fakeids, /save_cylinders
