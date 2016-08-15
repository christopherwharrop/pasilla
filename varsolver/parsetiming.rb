#!/usr/bin/ruby

class GPTLReader

  def initialize(fname)

    @total_time=0  # Total time of the program
    @records=[]    # Array of timing records in the timing data file
    
    # Loop over lines in the file to extract out the timing data
    data_start=false
    IO.readlines(fname).each { |line|

      # Process the line if we have seen the start of the timing data
      if data_start

        # Time to stop reading input if we hit a blank line
        break if line=~/^$/

        # Skip lines consisting of only blanks
        next if line=~/^\s+$/

        # Match the fields we want (indentation, routine, calls, wall time) 
        m=/^(\**\s+)(\S+)\s+(\S+)\s+\S+\s+(\S+)/.match(line)

        # Make a timing record for it
        record = {"indent" =>  m[1], "level" => m[1].length/2, "routine" => m[2], "calls" => m[3].to_i, "time" => m[4].to_f }

        # If this is the first line of data set the total time to the wall time of the record
        if @total_time == 0
          @total_time = record["time"]
        end

        # Add this record to our data array
        @records << record

      # Otherwise check to see if this line is the start of the timing data
      elsif line=~/Called  Recurse Wallclock max       min       self_OH  parent_OH/
        data_start=true
      end
    }

  end


  def print_report(level,threshold=0,routine=nil)

    # Get the records matching level and threshold constraints
    print_records = @records.find_all { |record| (record["level"] <= level) && (record["time"] / @total_time * 100.0 >= threshold) }

    # Extract records matching ancestry and descendants for routine name if specified
    unless routine.nil?

      # Collect the indices of all records matching the routine name
      routine_indices = []
      print_records.each_with_index { |r,i| routine_indices << i if r["routine"]==routine }

      # Mark records that call or are called by routines that call or are called by the specified routine name
      routine_indices.each { |routine_index|
        routine_level = print_records[routine_index]["level"]
        current_level = routine_level
        print_records[routine_index]["tag"] = true

        # Starting with the record matching routine, mark ancestral line back to main
        routine_index.downto(0) { |i|
          if print_records[i]["level"] < current_level
            print_records[i]["tag"] = true 
            current_level = print_records[i]["level"]
          end
        }

        # Starting with record matchin routine, mark direct descendants 
        (routine_index+1).upto(print_records.length-1) { |i|
          break if print_records[i]["level"] <= routine_level
          print_records[i]["tag"] = true 
        }
      }

      # Set the records to print to those marked as ancestors or descendants of routine calls
      print_records = print_records.find_all { |record| record["tag"] == true }
    end

    # Compute the width of the first field to output to dynamically size the output formatting
    routine_width = print_records.map { |record| record["indent"].length + record["routine"].length }.max

    # Print the header
    printf "%-#{routine_width}s %7s %8s %9s\n", "Routine","Calls","Walltime","% of Tot"

    # Print the matching records
    print_records.each { |record|
      printf "%-#{routine_width}s %7d %8.1f %8.1f%%\n", record["indent"]+record["routine"], record["calls"], record["time"], record["time"] / @total_time * 100.0
    }

  end


end

require 'optparse'
depth=999999
threshold=0
routine=nil
fielname=nil

# Get the command line options
opts = OptionParser.new
opts.banner = "Usage: parsetiming [options] filename"
opts.on("-d","--depth DEPTH",Integer) { |val| depth=val }
opts.on("-t","--threshold THRESHOLD",Integer) { |val| threshold=val }
opts.on("-r","--routine ROUTINE",String) { |val| routine=val }
opts.on("-h","--help") { puts opts; exit }
filename=opts.parse(ARGV).join

# Make sure a filename was provided
if filename.empty?
  puts opts
  exit
end

# Parse the file and print the report
r=GPTLReader.new(filename)
r.print_report(depth,threshold,routine)

