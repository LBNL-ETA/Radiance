#!/bin/csh -f
rm -f summer1.amb
rholo -n 24 -f summer1.hdk summer.hif 'RIF= insummer.rif AMB=summer1.amb QUA=Hi' DISK=15
