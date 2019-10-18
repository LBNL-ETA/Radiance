#!/bin/csh -fe
# RCSid $Id: testsuncal.csh,v 1.1 2019/10/18 23:42:52 greg Exp $
#
# Test solar calculations
#
# Produces:
# LAT LON SM TIME DAY MONTH SALT2 SAZI2
# all angles are degrees, as taken by sun2.cal
cnt 1000 | rcalc -f sun2.cal -e YEAR=2020 \
		-e 'LAT=160*rand(.359*recno-10)-80' \
		-e 'LON=360*rand(3.561*recno+16.6)-180' \
		-e 'SM=floor((LON+7.5)/15)*15' \
		-e 'TIME=14*rand(-7.5858*recno-71)+5' \
		-e 'MONTH=floor(1+11.99*rand(recno*.785+5.5))' \
		-e 'ndays(m):select(m,31,28,31,30,31,30,31,31,30,31,30,31)' \
		-e 'DAY=floor(1+.99*rand(recno*-71+9)*ndays(MONTH))' \
		-e '$1=LAT;$2=LON;$3=SM;$4=TIME;$5=DAY;$6=MONTH;$7=SALT;$8=SAZI' \
	> tests1.txt
# Produces:
# LAT LON SM TIME DAY MONTH SALT2 SAZI2 SALT SAZI SALT-SALT2 SAZI-SAZI2
rcalc -f sun.cal -e 'DEG:PI/180' -e 'RLAT=$1*DEG;RLON=$2*DEG;RSM=$3*DEG' \
	-e 'TIME=$4;DAY=$5;MONTH=$6' \
	-e '$1=$1;$2=$2;$3=$3;$4=$4;$5=$5;$6=$6;$7=$7;$8=$8' \
	-e '$9=SALT;$10=SAZI;$11=SALT-$7;$12=SAZI-$8' tests1.txt \
	> tests2.txt
