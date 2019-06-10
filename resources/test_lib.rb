#!/usr/bin/ruby

# Radiance test library

# Use with cmake build system deployed by NREL
# Run from [build]/resources
# usage: [ruby]test_lib[.rb] (runs all tests)
# usage: [ruby]test_lib[.rb] -n [n] -t ['testname', 'all']


$stderr.sync = true
require 'optparse'
require 'etc'
# default options
flag    = false

tst  = "all"


tproc = Etc.nprocessors
nproc = tproc / 2 # half cores by default, override with -n 
list    = ["x", "y", "z"]

clean = false

# parse arguments
file = __FILE__
ARGV.options do |opts|
  opts.on("-t", "--test=val", String)   { |val| tst = val }
  opts.on("-n", "--nproc=val", Integer)  { |val| nproc = val }

  opts.on("-c", "--clean") { |val| clean = val} 
  opts.parse!
end

# print banner
warn "Radiance Test Library"

### Cleanup ###

if clean == true

    puts "cleaning up"

    clean_list = []

    globs = [
        '../test/renders/*.oct',
        '../test/renders/*.amb',
        '../test/renders/*_ill.dat',
        '../test/renders/blinds_ill?.dat',
        '../test/renders/*_*.hdr',
        '../test/renders/*.unf',
        './test/gen/gen*.rad'
    ]

    globs.each do |glob|
        add_list = Dir.glob(glob)
        (clean_list << add_list).flatten!
    end

    #*.[cg]pm{,.leaf}

    add_files = [
 
        '../test/renders/inst_rad.txt',
        '../test/renders/combined.rad', 
        '../test/renders/fmirror.mtx',
        '../test/gen/replmarks.rad',
        '../test/util/test.mtx',
        '../test/cal/cnt.txt',
        '../test/cal/rcalc.txt',
        '../test/cal/total.txt',
        '../test/cal/histo.txt',
        '../test/cal/rlam.txt'

    ]
    (clean_list << add_files).flatten!

    clean_list.each do |rmfile| 
        if File.exist?(rmfile)
            File.delete(rmfile)
            puts "deleted #{rmfile}"
        end
    end
    exit

end

## END Cleanup ##

# print opts
warn "test: #{tst.inspect}"
warn "using #{nproc.inspect} cores (of #{tproc} total)"

test_in = tst

test_total = 0
@test_pass = 0
@test_fail = 0

all_tests = [
    'test_xform' , 'test_rad' , 'test_rtrace' , 'test_vwright' , 'test_rcollate' , 'test_dctimestep',
    'test_getinfo' , 'test_rmtxop' , 'test_oconv' , 'test_rfluxmtx' , 'test_dielectric_def' , 
    'test_dielectric_fish', 'test_replmarks', 'test_gensky', 'test_gendaymtx','test_genbox', 
    'test_genrev', 'test_genworm', 'test_genprism', 'test_gensurf', 'test_cnt', 'test_rcalc',
    'test_total', 'test_histo', 'test_rlam'
]


# report pass/fail status

def rcpf
    if $?.exitstatus == 0
        @test_pass +=1

        puts "result: PASS"
    else
        @test_fail +=1        
        puts "result: FAIL"

    end
end

# Number of processes to use on tests that run multi-core
NPROC = nproc

# Image reduction for comparisons
RDU_PFILT = "pfilt -1 -r 1 -x 128 -y 128 -pa 1"

# Image comparison command
IMG_CMP = "radcompare -v -rms 0.07 -max 1.5"

## tests from test/util

def test_vwright

    Dir.chdir('../test/util') do
        # generate test input
        cmd = "vwright -vf test.vf 3.5 > vwright.txt"
        system(cmd)

        # compare to reference
        cmd = "radcompare -v ref/vwright.txt vwright.txt"
        system(cmd)

        # report pass/fail
        rcpf

        # cleanup
        # File.delete("vwright.mtx")
    end

end

def make_test_mtx

    Dir.chdir('../test/util') do

        cmd = "rcollate -t ../renders/ref/rfmirror.mtx > test.mtx"
        system(cmd)

    end

end

def test_rcollate

    make_test_mtx if !File.exist?("../test/util/test.mtx")

    Dir.chdir('../test/util') do
        cmd = "radcompare -v ref/test.mtx test.mtx"
        system(cmd)
    end

    rcpf

    #File.delete("#{t}/test.mtx")

end

def test_dctimestep

    make_test_mtx if !File.exist?("../test/util/test.mtx")

    Dir.chdir('../test/util') do
        cmd = "gensky 3 21 10:15PST +s -g .3 -g 2.5 -a 36 -o 124 | genskyvec -m 1 -c .92 1.03 1.2 | dctimestep '!rmtxop -ff -t test.mtx' > dctimestep.mtx"
        system(cmd)
        cmd = "radcompare -v ref/dctimestep.mtx dctimestep.mtx"
        system(cmd)
    end

    rcpf

    #File.delete("#{t}/dctimestep.mtx")

end

def test_getinfo

    make_test_mtx if !File.exist?("../test/util/test.mtx")
    
    Dir.chdir('../test/util') do
        cmd = "getinfo -a Guppies 'Fredo the Frog' < test.mtx | getinfo > getinfo.txt"
        system(cmd)    
        cmd = "radcompare -v ref/getinfo.txt getinfo.txt"
        system(cmd)
    end

    rcpf
    
    #File.delete("#{t}/getinfo.txt")

end

def test_rmtxop

    make_test_mtx if !File.exist?("../test/util/test.mtx")

    Dir.chdir('../test/util') do
        cmd = "rmtxop -ff -c .3 .9 .2 test.mtx -c .7 .2 .3 -t test.mtx > rmtxop.mtx"
        puts("command: #{cmd}")
        system(cmd)
        cmd = "radcompare -v ref/rmtxop.mtx rmtxop.mtx"
        puts("command: #{cmd}")
        system(cmd)
    end
    
    rcpf
    #File.delete("#{t}/rmtxop.mtx")

end

### END test/util ###


### test/renders ###

def inst_oct
    Dir.chdir('../test/renders') do
        cmd = "rad -v 0 inst.rif"
        puts("command: #{cmd}")
        system(cmd)
    end
end

def test_xform
    combined_rad if !File.exist?("../test/renders/combined.rad")
    Dir.chdir('../test/renders') do
        cmd = "radcompare -v -max 0.04 ref/combined.rad combined.rad"
        puts("command: #{cmd}")
        system(cmd)
        rcpf
    end
end

def combined_rad
    Dir.chdir('../test/renders') do
        cmd = "xform -f combined_scene.rad | grep -v '^[   ]*#' > combined.rad"
        puts("command: #{cmd}")
        system(cmd)
    end
end

def test_rad
    Dir.chdir('../test/renders') do
        cmd = 'rad -n -s -e inst.rif > inst_rad.txt'
        system(cmd)
        cmd = 'radcompare -v ref/inst_rad.txt inst_rad.txt'
        puts("command: #{cmd}")
        system(cmd)
        rcpf
    end
end

def test_oconv
    inst_oct if !File.exist?("../test/renders/inst.oct")
    Dir.chdir('../test/renders') do
        cmd="radcompare -v ref/inst.oct inst.oct"
        system(cmd)
    end
    rcpf
end

def gen_rtmirror_fish_hdr
    Dir.chdir('../test/renders') do
        puts "generating input file: rtmirror_fish.hdr"
        cmd = "rad -v 0 mirror.rif OPT=mirror.opt"
        system(cmd)
        cmd = "vwrays -ff -vf fish.vf -x 2048 -y 2048 | rtrace -n #{NPROC} @mirror.opt -ffc -x 2048 -y 2048 mirror.oct | pfilt -1 -e +3 -r .6 -x /2 -y /2 > rtmirror_fish.hdr"
        #     rm -f mirror.opt
        puts("command: #{cmd}")
        system(cmd)
    end
end

def test_rtrace
    gen_rtmirror_fish_hdr if !File.exist?("../test/renders/rtmirror_fish.hdr")
    Dir.chdir('../test/renders') do
        cmd = "#{RDU_PFILT} rtmirror_fish.hdr | #{IMG_CMP} -h ref/mirror_fish.hdr -"
        puts("command: #{cmd}")
        system(cmd)
        rcpf
    end
end

def test_rfluxmtx
    gen_rfmirror_mtx if !File.exist?("../test/renders/rfmirror.mtx")
    Dir.chdir('../test/renders') do
        cmd = 'radcompare -v -max .4 -rms .05 ref/rfmirror.mtx rfmirror.mtx'
        puts("command: #{cmd}")
        system(cmd)
        rcpf
    end
end

def gen_rfmirror_mtx 
    Dir.chdir('../test/renders') do
        cmd = "rfluxmtx -n #{NPROC} -ab 2 -lw 1e-4 mirror.rad dummysky.rad basic.mat diorama_walls.rad closed_end.rad front_cap.rad glass_pane.rad antimatter_portal.rad > rfmirror.mtx"
        puts("command: #{cmd}")
        system(cmd)
    end
end

def gen_dielectric_oct
    Dir.chdir('../test/renders') do
        cmd = 'rad -v 0 dielectric.rif'
        puts("command: #{cmd}")
        system(cmd)
    end
end

def test_dielectric_def     #ref/dielectric_def.hdr dielectric_def.hdr
    gen_dielectric_def if !File.exist?("../test/renders/dielectric_def.hdr")
    gen_dielectric_def_ref if !File.exist?("../test/renders/ref/dielectric_def.hdr")
    Dir.chdir('../test/renders') do
        cmd = "#{RDU_PFILT} dielectric_def.hdr | #{IMG_CMP} ref/dielectric_def.hdr -"
        puts("command: #{cmd}")
        system(cmd)
        rcpf
    end
end

def gen_dielectric_def_ref
    gen_dielectric_def if !File.exist?('../test/renders/dielectric_def.hdr')
    Dir.chdir('../test/renders') do
        cmd = '#{RDU_PFILT} dielectric_def.hdr > ref/dielectric_def.hdr'
        puts("command: #{cmd}")
        system(cmd)
    end
end

def gen_dielectric_def #dielectric.oct
    gen_dielectric_oct if !File.exist?('../test/renders/dielectric.oct')
    Dir.chdir('../test/renders') do
        cmd = 'rad -v def dielectric.rif'
        puts("command: #{cmd}")
        system(cmd)
    end
end


### Reference and test for dielectric view  fish ###

def test_dielectric_fish   # ref/dielectric_fish.hdr dielectric_fish.hdr
    gen_dielectric_fish_hdr if !File.exist?('../test/renders/dielectric_fish.hdr')
    gen_ref_dielectric_fish_hdr if !File.exist?('../test/renders/ref/dielectric_fish.hdr')
    Dir.chdir('../test/renders') do
        cmd = "#{RDU_PFILT} dielectric_fish.hdr | #{IMG_CMP} ref/dielectric_fish.hdr -"
        puts("command: #{cmd}")
        system(cmd)
        rcpf
    end
end

def gen_ref_dielectric_fish_hdr #dielectric.oct
    gen_dielectric_fish_hdr if !File.exist?('../test/renders/dielectric_fish.hdr')
    Dir.chdir('../test/renders') do
        cmd = "#{RDU_PFILT} dielectric_fish.hdr > ref/dielectric_fish.hdr"
        puts("command: #{cmd}")
        system(cmd)
    end
end

def gen_dielectric_fish_hdr
    gen_dielectric_oct if !File.exist?('../test/renders/dielectric.oct')
    Dir.chdir('../test/renders') do
        cmd = 'rad -v  fish dielectric.rif'
        puts("command: #{cmd}")
        system(cmd)
    end
end

### End dielectric-fish  tests


### END test/renders ###


### test/gen ###

def test_replmarks
    Dir.chdir('../test/gen') do
        cmd = "replmarks -s 1 -x dummy.rad rmod markers.rad | grep -v '^[   ]*#' > replmarks.rad"
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
    #rm -f replmarks.rad
end


def test_gensky
    Dir.chdir('../test/gen') do
        cmd = 'gensky 11 15 14:21EST +s -g .25 -t 3.5 -a 40.7128 -o 74.006 > NYC11-15-14-21.rad'
        puts("command: #{cmd}")
        system(cmd)
        
        cmd = 'radcompare ref/NYC11-15-14-21.rad NYC11-15-14-21.rad'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
    #rm -f NYC11-15-14-21.rad
end


def test_gendaymtx
    Dir.chdir('../test/gen') do
        cmd = 'gendaymtx -r 90 -m 1 -g .3 .2 .1 -c .9 .9 1.2 test.wea > test.smx'
        puts("command: #{cmd}")
        system(cmd)
        
        cmd = 'radcompare ref/test.smx test.smx'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
    #rm -f test.smx
end

def test_genbox
    Dir.chdir('../test/gen') do
        cmd = 'genbox tmat tbox 1 2 3 -i -b .1 > genbox.rad'
        puts("command: #{cmd}")
        system(cmd)
        
        cmd = 'radcompare ref/genbox.rad genbox.rad'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
    #rm -f genbox.rad
end

def test_genrev
    Dir.chdir('../test/gen') do
        cmd = "genrev tmat trev 'sin(2*PI*t)' '2+cos(2*PI*t)' 16 -s > genrev.rad"
        puts("command: #{cmd}")
        system(cmd)
        
        cmd = 'radcompare ref/genrev.rad genrev.rad'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
    #rm -f genrev.rad
end

def test_genworm
    Dir.chdir('../test/gen') do
        cmd = "genworm tmat tworm '0' '5*sin(t)' '5*cos(t)' '.4-(.5-t)*(.5-t)' 8 > genworm.rad"
        puts("command: #{cmd}")
        system(cmd)
        
        cmd = 'radcompare ref/genworm.rad genworm.rad'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
    #rm -f genworm.rad
end

def test_genprism
    Dir.chdir('../test/gen') do
        cmd = "genprism tmat tprism 8 0 5 4 5 4 4 24.5 4 24.5 3 30 1.5 30 22 0 22 -l 0 0 -1.5 -r .2 > genprism.rad"
        puts("command: #{cmd}")
        system(cmd)
        
        cmd = 'radcompare ref/genprism.rad genprism.rad'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
    #rm -f genprism.rad
end

def test_gensurf
    Dir.chdir('../test/gen') do
        cmd = "gensurf tmat tsurf '15.5+x(theta(s),phi(t))' '10.5+y(theta(s),phi(t))' '30.75+z(theta(s),phi(t))' 4 15 -f basin.cal -s > gensurf.rad"
        puts("command: #{cmd}")
        system(cmd)
        
        cmd = 'radcompare ref/gensurf.rad gensurf.rad'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
    #rm -f gensurf.rad
end


### END test/gen ###


### test/cal ###

def test_cnt
    gen_cnt_txt if !File.exist?('../test/cal/cnt.txt')
    Dir.chdir('../test/cal') do
        cmd = 'radcompare ref/cnt.txt cnt.txt'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
end

def gen_cnt_txt
    Dir.chdir('../test/cal') do        
        cmd = 'cnt 5 3 2 > cnt.txt'
        puts("command: #{cmd}")
        system(cmd)
    end
end

def gen_rcalc_txt
    Dir.chdir('../test/cal') do        
        cmd = "rcalc -o 'Test $${v1} $$(s1) $${v2}' -e 'v1=$$1*$$2;v2=($$2-$$1)*exp($$3)' -s s1=HEY cnt.txt > rcalc.txt"
        puts("command: #{cmd}")
        system(cmd)
    end
end

def test_rcalc
    gen_rcalc_txt if !File.exist?('../test/cal/rcalc.txt')
    Dir.chdir('../test/cal') do
        cmd = 'radcompare ref/rcalc.txt rcalc.txt'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
end

def gen_total_txt
    gen_cnt_txt if !File.exist?('../test/cal/cnt.txt')
    Dir.chdir('../test/cal') do
        build = [
            'total cnt.txt > total.txt',
            'total -l cnt.txt >> total.txt',
            'total -u cnt.txt >> total.txt',
            'total -m cnt.txt >> total.txt',
            'total -s2.5 cnt.txt >> total.txt',
            'total -3 -r cnt.txt >> total.txt'
        ]
        build.each do |cmd|
            system(cmd)
        end
    end
end

def test_total
    gen_total_txt if !File.exist?('../test/cal/total.txt')
    Dir.chdir('../test/cal') do
        cmd = 'radcompare ref/total.txt total.txt'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
end

def gen_histo_txt
    gen_total_txt if !File.exist?('../test/cal/total.txt')
    Dir.chdir('../test/cal') do
        cmd = 'histo 0 60 5 < total.txt > histo.txt'
        puts("command: #{cmd}")
        system(cmd)

    end
end

def test_histo
    gen_histo_txt if !File.exist?('../test/cal/histo.txt')
    Dir.chdir('../test/cal') do
        cmd = 'radcompare ref/histo.txt histo.txt'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
end

def gen_rlam_txt
    gen_total_txt if !File.exist?('../test/cal/total.txt')
    gen_cnt_txt if !File.exist?('../test/cal/cnt.txt')
    gen_histo_txt if !File.exist?('../test/cal/histo.txt')

    Dir.chdir('../test/cal') do
        cmd = 'rlam -in 5 total.txt cnt.txt histo.txt > rlam.txt'
        puts("command: #{cmd}")
        system(cmd)

    end
end

def test_rlam
    gen_rlam_txt if !File.exist?('../test/cal/rlam.txt')
    Dir.chdir('../test/cal') do
        cmd = 'radcompare ref/rlam.txt rlam.txt'
        puts("command: #{cmd}")
        system(cmd)

        rcpf
    end
end


### END test/cal ###

### END test methods ###


# call the test already

if test_in == "all"
    all_tests.each do |test_in|
        test_total += 1
        puts "running test: #{test_in}"
        method(test_in).call
        puts ''
    end
else
    test_total += 1
    puts "running test: #{test_in}"
    method(test_in).call
    puts ''
end

puts "### Total tests: #{test_total} (Passed: #{@test_pass} Failed: #{@test_fail}) ###"

