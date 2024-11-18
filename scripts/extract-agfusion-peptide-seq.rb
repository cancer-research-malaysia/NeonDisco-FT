#!/usr/bin/env ruby

# Check if the correct number of arguments is provided
if ARGV.length != 2
	puts "Usage: ruby extract-agfusion-peptide-seq.rb <input_directory> <output_directory>"
	exit 1
end
  
# Function to find all .fa files recursively
def find_fa_files(dir, files_array)
	Dir.foreach(dir) do |entry|
	  next if entry == '.' || entry == '..'
	  
	  full_path = File.join(dir, entry)
	  if File.directory?(full_path)
		find_fa_files(full_path, files_array)
	  elsif entry.end_with?('.fa')
		files_array << full_path
	  end
	end
end
  
def extract_identifier(filename)
	# Extract the identifier from the filename (GENEIDS_COORDINATES_sNUMBER)
	if filename =~ /(.+?)_protein\.fa$/
	  return $1
	else
	  return "unknown"
	end
end
  
def process_out_of_frame(sequence)
	asterisk_pos = sequence.index('*')
	return nil unless asterisk_pos
	
	start_pos = [asterisk_pos - 13, 0].max
	before_part = sequence[start_pos...asterisk_pos]
	
	sequence[start_pos..]
end
  
def process_in_frame(sequence)
	asterisk_pos = sequence.index('*')
	return nil unless asterisk_pos
	
	start_pos = [asterisk_pos - 13, 0].max
	before_part = sequence[start_pos...asterisk_pos]
	
	after_part = sequence[(asterisk_pos + 1)...(asterisk_pos + 14)]
	after_part = after_part || ""
	
	before_part + '*' + after_part
end
  
def process_fa_file(file_path)
	out_of_frame_peptides = []
	in_frame_peptides = []
	current_header = nil
	current_sequence = ""
	
	File.readlines(file_path).each do |line|
	  line.chomp!
	  
	  if line.start_with?(">")
		  if current_header && !current_sequence.empty?
		    if current_header.include?("effect: out-of-frame")
			    processed_sequence = process_out_of_frame(current_sequence)
			    out_of_frame_peptides << processed_sequence if processed_sequence
		    elsif current_header.include?("effect: in-frame")
			    processed_sequence = process_in_frame(current_sequence)
			    in_frame_peptides << processed_sequence if processed_sequence
		    end
		    current_sequence = ""
		  end
		  current_header = line
	  else
		  current_sequence += line
	  end
	end
	
	# Process the last sequence
	if current_header && !current_sequence.empty?
	  if current_header.include?("effect: out-of-frame")
		  processed_sequence = process_out_of_frame(current_sequence)
		  out_of_frame_peptides << processed_sequence if processed_sequence
	  elsif current_header.include?("effect: in-frame")
		  processed_sequence = process_in_frame(current_sequence)
		  in_frame_peptides << processed_sequence if processed_sequence
	  end
	end
	
	# Get unique sequences only
  return out_of_frame_peptides.uniq, in_frame_peptides.uniq
end
  
# Get base directory from command line argument
base_dir = ARGV[0]
out_dir = ARGV[1]

# Initialize arrays for collecting fa files
fa_files = []

# Find all .fa files
find_fa_files(base_dir, fa_files)

# Print total number of files found
puts "Found #{fa_files.length} .fa files to process..."
# puts "#{fa_files}"
puts "Processing files..."

# Initialize arrays for cumulative unique peptides
all_out_of_frame = []
all_in_frame = []

# Process each file and save results separately
total_out_of_frame = 0
total_in_frame = 0

fa_files.each_with_index do |file_path, index|
  puts "Processing file #{index + 1}/#{fa_files.length}: #{file_path}"

  # Extract identifier from filename
  filename = File.basename(file_path)
  identifier = extract_identifier(filename)

  # Process the file
  out_of_frame, in_frame = process_fa_file(file_path)

  # Update totals
  total_out_of_frame += out_of_frame.length
  total_in_frame += in_frame.length

  # Add to cumulative arrays
  all_out_of_frame.concat(out_of_frame)
  all_in_frame.concat(in_frame)

  # Write results to separate files for this input file
  File.open(File.join(out_dir, "#{identifier}_out_of_frame_peptides.txt"), 'w') do |file|
    out_of_frame.each do |sequence|
    file.puts sequence
    end
  end

  File.open(File.join(out_dir, "#{identifier}_in_frame_peptides.txt"), 'w') do |file|
    in_frame.each do |sequence|
    file.puts sequence
    end
  end

  # Print individual file statistics
  puts "  Found #{out_of_frame.length} unique out-of-frame peptides"
  puts "  Found #{in_frame.length} unique in-frame peptides"
end

# Get unique cumulative peptides
all_out_of_frame.uniq!
all_in_frame.uniq!

# Write cumulative results
File.open(File.join(out_dir, "ALL_out_of_frame_unique.txt"), 'w') do |file|
  all_out_of_frame.each do |sequence|
    file.puts sequence
  end
end

File.open(File.join(out_dir, "ALL_in_frame_unique.txt"), 'w') do |file|
  all_in_frame.each do |sequence|
    file.puts sequence
  end
end

# Print summary
puts "\nProcessing complete!"
puts "Individual files:"
puts "  Found #{total_out_of_frame} total out-of-frame peptides"
puts "  Found #{total_in_frame} total in-frame peptides"
puts "\nCumulative unique peptides:"
puts "  Found #{all_out_of_frame.length} unique out-of-frame peptides across all files"
puts "  Found #{all_in_frame.length} unique in-frame peptides across all files"
puts "Results saved in: #{out_dir}"
