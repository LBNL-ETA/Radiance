#!/usr/local/bin/python

import sys

cal_units = [
	'test_calc',
	'test_cnt',
	'test_ev',
	'test_histo',
	'test_lam',
	'test_neat',
	'test_rcalc',
	'test_tabfunc',
	'test_total',
]

cv_units = [
	'test_3ds2mgf',
	'test_arch2rad',
	'test_ies2rad',
	'test_lampcolor',
	'test_mgf2inv',
	'test_mgf2rad',
	'test_nff2rad',
	'test_obj2rad',
	'test_thf2rad',
	'test_tmesh2rad',
	'test_rad2mgf',
	'test_mgf2meta',
	'test_mgfilt',
]

gen_units = [
	'test_genbeads',
	'test_genblinds',
	'test_genbox',
	'test_genbranch',
	'test_gencat',
	'test_genclock',
	'test_genmarble',
	'test_genprism',
	'test_genrev',
	'test_gensky',
	'test_gensurf',
	'test_genworm',
	'test_replmarks',
	'test_mkillum',
	'test_xform',
]

hd_units = [
	'test_rhcopy',
	'test_rhinfo',
	'test_rholo',
	'test_rhpict',
	#'test_genrhenv',
	'test_genrhgrid',
	# XXX device drivers?
]

meta_units = [
	'test_cv',
	'test_meta2tga',
	'test_pexpand',
	'test_plot4',
	'test_plotin',
	'test_psort',
	'test_psmeta',
	'test_gcomp',
	'test_bgraph',
	'test_dgraph',
	'test_igraph',
	#'test_x11meta',
]

meta_special_units = [
	'test_mt160',
	'test_mt160l',
	'test_mtext',
	'test_okimate',
	'test_mx80',
	'test_imagew',
	'test_impress',
	'test_aed5',
	'test_tcurve',
	'test_tbar',
	'test_tscat',
	'test_plotout',
]

ot_units = [
	'test_oconv',
	'test_getbbox',
	'test_obj2mesh',
]

px_units = [
	'test_macbethcal',
	'test_normtiff',
	'test_oki20',
	'test_oki20c',
	'test_pcomb',
	'test_pcompos',
	'test_pcond',
	'test_pcwarp',
	'test_pextrem',
	'test_pfilt',
	'test_pflip',
	'test_pinterp',
	'test_protate',
	'test_psign',
	'test_pvalue',
	'test_ra_avs',
	'test_ra_bn',
	'test_ra_gif',
	'test_ra_hexbit',
	'test_ra_pict',
	'test_ra_ppm',
	'test_ra_pr',
	'test_ra_pr24',
	'test_ra_ps',
	'test_ra_rgbe',
	'test_ra_t16',
	'test_ra_t8',
	'test_ra_xyze',
	'test_ra_tiff',
	'test_t4027',
	'test_ttyimage',
	#'test_ximage',
	#'test_xshowtrace',
]

px_special_units = [
	'test_ra_im',
	'test_psum',
	'test_t4014',
	'test_paintjet',
	'test_mt160r',
	'test_greyscale',
	'test_colorscale',
	'test_d48c',
]

rt_units = [
	'test_lookamb',
	'test_rpict',
	'test_rtrace',
	#'test_rview',
]

util_units = [
	'test_rad',
	'test_findglare',
	'test_glarendx',
	'test_rpiece',
	'test_vwrays',
	'test_vwright',
	'test_getinfo',
	'test_makedist',
	#'test_xglaresrc',
	'test_glrad',
	'test_ranimate',
	'test_ranimove',
]


def run_tests(unitgroup):
	print '---- unit group %s ----' % unitgroup
	for unit in globals()[unitgroup + '_units']:
		try:
			mod = __import__('py_tests.'+unit,globals(),locals(),['py_tests'])
			print '%-18s' % unit,
			sys.stdout.flush()
			mod.main()
		except ImportError, e:
			#raise
			pass

def main():
	run_tests('cal')
	run_tests('cv')
	run_tests('gen')
	run_tests('hd')
	run_tests('meta')
	run_tests('ot')
	run_tests('px')
	run_tests('rt')
	run_tests('util')

if __name__ == '__main__':
	main()

