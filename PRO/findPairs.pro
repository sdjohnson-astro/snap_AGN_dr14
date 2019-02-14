PRO findPairs, quasars, pairs
	
	quasars = mrdfits('../catalogs/DR14Q_v4_4_Paris2017.fits', 1)
	
	; remove bad quasars
	index = where(quasars.plate ne 7306 AND quasars.MJD ne 56684 AND quasars.fiberID ne 938)
	quasars = quasars[index]
	
	quasars = quasars[where(quasars.z gt 0.36)]
	fibers = mrdfits('../catalogs/DR14Q_spec1_sj211.fits', 1)
	
	struct_add_field, quasars, 'm_FUV', 0d
	struct_add_field, quasars, 'm_NUV', 0d
	struct_add_field, quasars, 'm_u', 0d
	struct_add_field, quasars, 'm_g', 0d
	struct_add_field, quasars, 'm_r', 0d
	struct_add_field, quasars, 'm_i', 0d
	struct_add_field, quasars, 'm_z', 0d
	struct_add_field, quasars, 'SNMEDIAN_u', 0d
	struct_add_field, quasars, 'SNMEDIAN_g', 0d
	struct_add_field, quasars, 'SNMEDIAN_r', 0d
	struct_add_field, quasars, 'SNMEDIAN_i', 0d
	struct_add_field, quasars, 'SNMEDIAN_z', 0d
	
	quasars.m_u = reform(quasars.psfMag[0, *])
	quasars.m_g = reform(quasars.psfMag[1, *])
	quasars.m_r = reform(quasars.psfMag[2, *])
	quasars.m_i = reform(quasars.psfMag[3, *])
	quasars.m_z = reform(quasars.psfMag[4, *])
	
	quasars.m_FUV = 22.5 - 2.5*alog10(quasars.FUV);-2.5*alog10(quasars.FUV*1e-6/3631d)
	quasars.m_NUV = 22.5 - 2.5*alog10(quasars.NUV);-2.5*alog10(quasars.NUV*1e-6/3631d)
	
	
	; Match the two tables
	spherematch, quasars.ra, quasars.dec, fibers.ra, fibers.dec, 1d/3600d, match1, match2
	
	quasars = quasars[match1]
	fibers = fibers[match2]
	
	quasars.SNMEDIAN_u = fibers.SNMEDIAN_u
	quasars.SNMEDIAN_g = fibers.SNMEDIAN_g
	quasars.SNMEDIAN_r = fibers.SNMEDIAN_r
	quasars.SNMEDIAN_i = fibers.SNMEDIAN_i
	quasars.SNMEDIAN_z = fibers.SNMEDIAN_z
	
   mwrfits, quasars, '../catalogs/quasasrs.fits', /create
	
	print, 'Catalog complete'
	
	pair = {index:0l, name_fg:'', ra_fg:0d, dec_fg:0d, plate_fg:0l, MJD_fg:0l, fiberID_fg:0l, redshift_fg:0d, $
		     m_fuv_fg:0d, m_nuv_fg:0d, m_u_fg:0d, m_g_fg:0d, m_r_fg:0d, m_i_fg:0d, m_z_fg:0d, $
			  SNMEDIAN_u_fg:0d, SNMEDIAN_g_fg:0d, SNMEDIAN_r_fg:0d, SNMEDIAN_i_fg:0d, SNMEDIAN_z_fg:0d, $
			  name_bg:'', ra_bg:0d, dec_bg:0d, plate_bg:0l, MJD_bg:0l, fiberID_bg:0l, redshift_bg:0d, $
			  m_fuv_bg:0d, m_nuv_bg:0d, m_u_bg:0d, m_g_bg:0d, m_r_bg:0d, m_i_bg:0d, m_z_bg:0d, $
			  SNMEDIAN_u_bg:0d, SNMEDIAN_g_bg:0d, SNMEDIAN_r_bg:0d, SNMEDIAN_i_bg:0d, SNMEDIAN_z_bg:0d, $
			  wave_Lya_forest:0d, wave_CII:0d, wave_CIV:0d, wave_FeII:0d, wave_MgII:0d, $
			  quality_fg:-1l, quality_bg:-1l, $
			  dv:0d, theta:0d, d:0d}
	
	
	print, 'Finding pairs...'
	
	spherematch, quasars.ra, quasars.dec, quasars.ra, quasars.dec, 5d/60d, match1, match2, maxmatch=0
	
	quasars1 = quasars[match1]
	quasars2 = quasars[match2]
	
	pairs = replicate(pair, n_elements(match1))
	pairs.name_fg = quasars1.SDSS_NAME
	pairs.ra_fg = quasars1.ra
	pairs.dec_fg = quasars1.dec
	pairs.plate_fg = quasars1.plate
	pairs.MJD_fg = quasars1.MJD
	pairs.fiberID_fg = quasars1.fiberID
	pairs.redshift_fg = quasars1.z
	pairs.m_fuv_fg = quasars1.m_fuv
	pairs.m_nuv_fg = quasars1.m_nuv
	pairs.m_u_fg = quasars1.m_u
	pairs.m_g_fg = quasars1.m_g
	pairs.m_r_fg = quasars1.m_r
	pairs.m_i_fg = quasars1.m_i
	pairs.m_z_fg = quasars1.m_z
	
	pairs.SNMEDIAN_u_fg = quasars1.SNMEDIAN_u
	pairs.SNMEDIAN_g_fg = quasars1.SNMEDIAN_g
	pairs.SNMEDIAN_r_fg = quasars1.SNMEDIAN_r
	pairs.SNMEDIAN_i_fg = quasars1.SNMEDIAN_i
	pairs.SNMEDIAN_z_fg = quasars1.SNMEDIAN_z
	
	pairs.name_bg = quasars2.SDSS_NAME
	pairs.ra_bg = quasars2.ra
	pairs.dec_bg = quasars2.dec
	pairs.plate_bg = quasars2.plate
	pairs.MJD_bg = quasars2.MJD
	pairs.fiberID_bg = quasars2.fiberID
	pairs.redshift_bg = quasars2.z
	pairs.m_fuv_bg = quasars2.m_fuv
	pairs.m_nuv_bg = quasars2.m_nuv
	pairs.m_u_bg = quasars2.m_u
	pairs.m_g_bg = quasars2.m_g
	pairs.m_r_bg = quasars2.m_r
	pairs.m_i_bg = quasars2.m_i
	pairs.m_z_bg = quasars2.m_z
	
	pairs.SNMEDIAN_u_bg = quasars2.SNMEDIAN_u
	pairs.SNMEDIAN_g_bg = quasars2.SNMEDIAN_g
	pairs.SNMEDIAN_r_bg = quasars2.SNMEDIAN_r
	pairs.SNMEDIAN_i_bg = quasars2.SNMEDIAN_i
	pairs.SNMEDIAN_z_bg = quasars2.SNMEDIAN_z
	
	pairs.wave_Lya_forest = 1215.67*(1 + pairs.redshift_bg)
	pairs.wave_CII = 1334.53*(1 + pairs.redshift_fg)
	pairs.wave_CIV = 1548.19*(1 + pairs.redshift_fg)
	pairs.wave_FeII = 2600.17*(1 + pairs.redshift_fg)
	pairs.wave_MgII = 2796.35*(1 + pairs.redshift_fg)
	pairs.dv = calcvdiff(pairs.redshift_bg, pairs.redshift_fg)
	pairs.theta = coordtotheta(pairs.ra_fg, pairs.dec_fg, pairs.ra_bg, pairs.dec_bg, /d)
	pairs.d = thetatorho(pairs.theta, pairs.redshift_fg)
	
	
	index = where(pairs.dv lt -10000d AND pairs.theta gt 1d AND pairs.d lt 500d)
	pairs = pairs[index]
	
	index = where(pairs.wave_lya_forest lt pairs.wave_MgII)
	pairs = pairs[index]
	
	pairs = pairs[sort(pairs.d)]
	
	pairs.index = lindgen(n_elements(pairs))
	
	
	mwrfits, pairs, '../catalogs/pairs.fits', /create
	
	
	
	
	
END