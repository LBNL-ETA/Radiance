#!/bin/csh -f
goto cont
oconv -i summercabin.oct -f ext1.rad > ext1sum.oct
rad -v 'cor -vf vf/corner' -w insummer.rif PIC=ext1sum OCT=ext1sum.oct EXP=-1
rm ext1sum.oct
oconv -i summercabin.oct -f ext2.rad > ext2sum.oct
rad -v 'cor -vf vf/corner' -w insummer.rif PIC=ext2sum OCT=ext2sum.oct EXP=-1
rm ext2sum.oct
cont:
oconv -f summer0_env.rad ext1.rad > ext1env.oct
rad -v 'cor -vf vf/corner' -w insummer.rif QUA=Hi PIC=ext1env \
	RAW=ext1raw ZF=ext1raw OCT=ext1env.oct EXP=-1 'render=-ab 1'
rm ext1env.oct
oconv -f summer0_env.rad ext2.rad > ext2env.oct
rad -v 'cor -vf vf/corner' -w insummer.rif QUA=Hi PIC=ext2env \
	RAW=ext2raw ZF=ext2raw OCT=ext2env.oct EXP=-1 'render=-ab 1'
rm ext2env.oct
