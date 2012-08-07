#!/usr/bin/env ruby
# RCSid: $Id: objview.rb,v 2.1 2012/08/07 18:17:17 greg Exp $
#
# objview.rb, v1.1 2012.07.04 RPG
# Ruby port of objview.csh
#
# Make a nice view of an object, or luminaire
# Arguments are scene input files
#

require 'optparse'
require 'ostruct'
require 'tempfile'

class Optparse
  #
  # Return a structure describing the options.
  #
  def self.parse(args)
    # The options specified on the command line will be collected in *options*.
    options = OpenStruct.new
 
    opts = OptionParser.new do |opts|
      opts.banner = "Usage: (ruby) objview.rb [options] [radiance scene input files]"

      opts.separator ""
      opts.separator "Options:"
      
      options.radopt = []
      options.glradopt = []
      
      options.up = "Z"
      opts.on("-u, --up", String, [:X, :Y, :Z], "View up vector (X, Y, or Z)") do |n|
        options.up = n
      end
      
      #Specify display device (qt is only option for Windows)
      if /mswin/.match(RUBY_PLATFORM) or /mingw/.match(RUBY_PLATFORM)
        options.dev = "qt"
      else
        options.dev = "x11"
      end
      opts.on("-o, --odev", String, [:x11, :qt], "Output device (x11 or qt)") do |n|
        options.dev = n
        options.radopt << "-o #{options.dev}"
      end
      
      options.procs = "1"
      opts.on("-N, --nprocs", String, "Number of processors (UNIX ony)") do |n|
        options.procs = n
        options.radopt << "-N #{options.procs}"
      end
      
      options.vw = "XYZ"
      opts.on("-v", "--view", String, "View specification") do |n|
        options.view = n
        options.radopt << "-v #{options.view}"
      end
     
      options.ltv = false
      opts.on("-l", "--ltview", "Luminaire viewer mode") do |n|
        options.ltv = true
        options.vw = "y"
      end
      
      #set "room" size for ltview, default is 4' square box (in meters)
      options.rsiz = "0.6096"
      opts.on("-r", "--roomsize", String, "Room scaling factor for ltview") do |n|
        options.rsiz = n
      end
      
      opts.on("-d", "--debug", "Run verbosely") do |n|
        options.verbose = true
      end
      
      opts.on("-g", "--ogl", "Use OpenGL (glrad)") do |n|
        options.ogl = true
      end
      
      opts.on("-e", "--rad", "Set rad options") do |n|
        options.radopt << n
      end
      
      opts.on("-S", "--glrad", "Set glrad options") do |n|
        options.glradopt = "-S #{n}"
      end     
      
      # No argument, shows at tail.  This will print an options summary.
      opts.on_tail("-h", "--help", "Show this message") do
        puts opts
        exit
      end
      
    end
    
    opts.parse!(args)
    options
    
  end  # parse()
  
end  # class Optparse

options = Optparse.parse(ARGV)

# utility function: print statement and execute as system call
def exec_statement(s)
  if /mswin/.match(RUBY_PLATFORM) or /mingw/.match(RUBY_PLATFORM)
    s = s.gsub("/", "\\")
  end
  puts "'#{s}'"
  result = system(s)
  puts
  return result
end

exec = "system"
if options.verbose == true
  puts options
  puts "Input file(s): #{ARGV}"
  exec = "exec_statement"
end

#not needed(?)
#octree = Tempfile.new(['objview_', '.oct'])
#ambf = Tempfile.new(['objview_', '.amb'])

if options.ltv == true
  rendopts = "-ab 1 -ds .15"
end

# not sure how these options are supposed to work... RPG
#
# 	case -s:
# 	case -w:
# 		set opts=($opts $argv[1])
# 		breaksw
# 	case -b*:
# 		set rendopts=($rendopts -bv)
# 		breaksw
# 	case -V:
#


if not ARGV[0]
  puts "No input files specified"
  exit 1
end
  
if options.ogl == true
  if rendopts
    puts "bad option for glrad"
    exit 1
  else
    if options.glradopt == true 
      puts "bad option for rad"
      exit 1
    end
  end  
end

if options.ltv == false
  # create lights input file for objview
  lights = Tempfile.new(['lt_', '.rad'])
  lights.write("void glow dim 0 0 4 .1 .1 .15 0\n")
  lights.write("dim source background 0 0 4 0 0 1 360\n")
  lights.write("void light bright 0 0 3 1000 1000 1000\n")
  lights.write("bright source sun1 0 0 4 1 .2 1 5\n")
  lights.write("bright source sun2 0 0 4 .3 1 1 5\n")
  lights.write("bright source sun3 0 0 4 -1 -.7 1 5\n")
  if options.verbose == true
    lights.rewind
    lightsread = lights.read
    puts "\nlights file:"
    puts lightsread
  end
  lights.close
end

if options.ltv == true
  #create room input file for ltview
  room = Tempfile.new(['lt_', '.rad'])
  room.write("void plastic surf\n0\n0\n5\n.2 .2 .2 0 0\n\n")
  room.write("surf polygon floor\n0\n0\n12\n")
  room.write("-#{options.rsiz} #{options.rsiz} -#{options.rsiz}\n")
  room.write("-#{options.rsiz} -#{options.rsiz} -#{options.rsiz}\n")
  room.write("#{options.rsiz} -#{options.rsiz} -#{options.rsiz}\n")
  room.write("#{options.rsiz} #{options.rsiz} -#{options.rsiz}\n\n")

  room.write("surf polygon clg\n0\n0\n12\n")
  room.write("-#{options.rsiz} #{options.rsiz} #{options.rsiz}\n")
  room.write("-#{options.rsiz} -#{options.rsiz} #{options.rsiz}\n")
  room.write("#{options.rsiz} -#{options.rsiz} #{options.rsiz}\n")
  room.write("#{options.rsiz} #{options.rsiz} #{options.rsiz}\n\n")

  room.write("surf polygon wall-north\n0\n0\n12\n")
  room.write("-#{options.rsiz} #{options.rsiz} -#{options.rsiz}\n")
  room.write("#{options.rsiz} #{options.rsiz} -#{options.rsiz}\n")
  room.write("#{options.rsiz} #{options.rsiz} #{options.rsiz}\n")
  room.write("-#{options.rsiz} #{options.rsiz} #{options.rsiz}\n\n")

  room.write("surf polygon wall-east\n0\n0\n12\n")
  room.write("#{options.rsiz} #{options.rsiz} -#{options.rsiz}\n")
  room.write("#{options.rsiz} -#{options.rsiz} -#{options.rsiz}\n")
  room.write("#{options.rsiz} -#{options.rsiz} #{options.rsiz}\n")
  room.write("#{options.rsiz} #{options.rsiz} #{options.rsiz}\n\n")  
  
  room.write("surf polygon wall-west\n0\n0\n12\n")
  room.write("-#{options.rsiz} #{options.rsiz} -#{options.rsiz}\n")
  room.write("-#{options.rsiz} #{options.rsiz} #{options.rsiz}\n")
  room.write("-#{options.rsiz} -#{options.rsiz} #{options.rsiz}\n")
  room.write("-#{options.rsiz} -#{options.rsiz} -#{options.rsiz}\n\n")
  if options.verbose == true
    room.rewind
    roomread = room.read
    puts "\nroom file:"
    puts roomread
  end
  room.close
end

# create .rif file
rif = Tempfile.new(['ov_', '.rif'])
if options.ltv == true
  rif.write("scene= #{ARGV.join(' ')} #{room.path}\n")
else
  rif.write("scene= #{ARGV.join(' ')} #{lights.path}\n")
end
rif.write("EXPOSURE= .5\n")
rif.write("UP= #{options.up}\n")
rif.write("view= #{options.vw}\n")
rif.write("oconv= -f\n")
rif.write("render= #{rendopts}\n")
#doesn't rad handle the construction and destruction of these files (below) on its own? RPG
#rif.write("OCTREE= #{octree.path}\n")
#rif.write("AMBF= #{ambf.path}\n")
if options.verbose == true
  rif.rewind
  rifread = rif.read
  puts "\nrif file:"
  puts rifread
end
rif.close

puts "radopt = #{options.radopt.join(" ")}"

#run executive control program of choice
if options.ogl == true
  exec("glrad #{options.radopt.join(" ")} #{rif.path}")
else
  exec("rad -o #{options.dev} #{options.radopt.join(" ")} #{rif.path}")
end

# clean up temp files
if coptions.ltv == true
  room.unlink
else  
  lights.unlink
end
rif.unlink
#ambf.unlink
#octree.unlink

exit 0

