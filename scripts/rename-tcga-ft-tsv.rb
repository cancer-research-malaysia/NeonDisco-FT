#!/usr/bin/env ruby

require 'fileutils'

# Check if the correct number of arguments is provided
if ARGV.length != 4
  puts "Usage: ruby wrangle-tcga-ft-tsv.rb <base_directory_path> <subdirectory> <uuid_text_file_path> <output_directory>"
  exit 1
end

###### Define functions #######
def find_tsv_file(tool_dir_path, id, tsv_file_name, dir_suffix)
  dir_path = File.join(tool_dir_path, "#{id}#{dir_suffix}")
  file_path = File.join(dir_path, tsv_file_name)
  File.exist?(file_path) ? file_path : nil
end

def read_ids(file_path)
  File.readlines(file_path).map(&:strip)
end

def copy_and_rename_file(sample_tsv_filepath, id, tool, out_dir)
  case tool
  when 'Arriba'
    new_name = "#{id}_arr.tsv"
    tool_dir = File.join(out_dir, tool)
  when 'FusionCatcher'
    new_name = "#{id}_fc.tsv"
    tool_dir = File.join(out_dir, tool)
  else
    raise ArgumentError, "Unsupported tool: #{tool}"
  end

  # Create the tool-specific directory
  FileUtils.mkdir_p(tool_dir)
  ##new_path = File.join(File.dirname(File.dirname(sample_tsv_filepath)), new_name)
  out_path = File.join(tool_dir, new_name)
  FileUtils.cp(sample_tsv_filepath, out_path)
  puts "Copied and renamed: #{out_path}"
    
  rescue StandardError => e
    puts "Error processing #{id}: #{e.message}"
end

############## Main execution ##############
# Configuration from command-line arguments
base_dir = ARGV[0]
tool_subdirectory = ARGV[1]
uuid_file_path = ARGV[2]
out_dir = ARGV[3]

# check what is base_dir value
if tool_subdirectory == 'Arriba'
    tsv_file_name = 'arriba-fusions.tsv'
    dir_suffix = '_arr_outdir'
elsif tool_subdirectory == 'FusionCatcher'
    tsv_file_name = 'final-list_candidate-fusion-genes.txt'
    dir_suffix = '_fc_outdir'
else
    puts "Error: Unsupported tool subdirectory. Must be either 'Arriba' or 'FusionCatcher'."
    exit 1
end

# Validate inputs
input_dir_path = File.join(base_dir, tool_subdirectory)

unless File.directory?(input_dir_path)
  puts "Error: Invalid tool subdirectory or base directory path."
  exit 1
end

unless File.exist?(uuid_file_path)
  puts "Error: UUID text file not found."
  exit 1
end

puts "Base directory: #{base_dir}"
puts "Tool subdirectory: #{tool_subdirectory}"
puts "UUID file: #{uuid_file_path}"

ids = read_ids(uuid_file_path)
ids.each do |id|
  tsv_file_path = find_tsv_file(input_dir_path, id, tsv_file_name, dir_suffix)
  if tsv_file_path
    copy_and_rename_file(tsv_file_path, id, tool_subdirectory, out_dir)
  else
    puts "TSV file not found for ID: #{id}"
  end
end