#!/usr/bin/env ruby

framework_root_directory = File.absolute_path(File.dirname __FILE__)
$LOAD_PATH << "#{framework_root_directory}"

require "fileutils"
require "framework/xmlparser"
require "framework/external-script"
require "framework/environment"
require "framework/logger"
require "shellwords"

def generateCCPLConfiguration(case_name, models)
  models.each do |model|
    FileUtils.mkdir_p("CCPL_dir/config/all")
    model_name = model["name"]
    running_settings = model["overridable_settings"]["running_settings"]
    file_path = File.absolute_path("CCPL_dir/config/all/env_run.xml")
    open(file_path, "w") do |f|
      f.puts "<?xml version=\"1.0\" ?>"
      f.puts "<Time_setting"
      f.puts "  case_name=\"#{case_name || "unknown"}\""
      f.puts "  model_name=\"#{case_name || "unknown"}\""
      f.puts "  run_type=\"#{running_settings["run.type"] || "initial"}\""
      f.puts "  leap_year=\"#{running_settings["run"]["leap_year"] || "on"}\""
      f.puts "  start_date=\"#{running_settings["run"]["start_date"] || "00000101"}\""
      f.puts "  start_second=\"#{running_settings["run"]["start_second"] || "0"}\""
      f.puts "  reference_date=\"#{running_settings["run"]["reference_date"] || "00000101"}\""
      f.puts "  rest_freq_unit=\"#{running_settings["write_restart"]["freq_unit"] || "none"}\""
      f.puts "  rest_freq_count=\"#{running_settings["write_restart"]["freq_count"] || "1"}\""
      f.puts "  stop_option=\"#{running_settings["run"]["stop_option"] || "day"}\""
      f.puts "  stop_date=\"#{running_settings["run"]["stop_date"] || "00000101"}\""
      f.puts "  stop_second=\"#{running_settings["run"]["stop_second"] || "0"}\""
      f.puts "  stop_n=\"#{running_settings["run"]["stop_n"] || "-999"}\""
      f.puts "/>"
    end
  end
end

ScriptLogger.getLogDirectory(false)

# Parse paramters

# Load configuration
case_configuration = XmlParser::CaseConfiguration.new("config/case.xml")
models = case_configuration.getModelsFlatList

models.each() do |model|
  path = File.absolute_path("run/" + model["type"] + "/" + model["name"])
  path << "/#{model["ensemble_id"]}" if model["ensemble_id"] != nil
  model[".run_path"] = path
end

# Set some basic environment
basic_env = Hash.new
basic_env["CASE_NAME"] = case_configuration.case_name

case_configuration.machine_options.each do |key, value|
  basic_env["MACHINE_OPTIONS_#{key.upcase}"]  = value
end if case_configuration.machine_options != nil

generateCCPLConfiguration(case_configuration.case_name, models)

tempfile = Tempfile.new("submit")
models.each do |model|
  pre_submit_script = File.absolute_path(
    "config/models/#{model["type"]}/#{model["name"]}/pre_submit.sh")
  new_env, = EnvironmentSetting.setEnvironment(basic_env, model)
  status, stdout, = ExternalScript.callScript(pre_submit_script, new_env, true) if File.exist?(pre_submit_script)
  submitSettings = model["overridable_settings"]["submit_settings"]
  model_cmd = Array.new
  model_cmd.push("#{model[".run_path"]}/exe/#{model["name"]}")
  if submitSettings != nil then
    submitSettings.each do |paramter|
      param_str = paramter
      param_str.gsub!(/\$\{ENSEMBLE_ID\}/, model["ensemble_id"].to_s)
      param_str.gsub!(/\$\{ENSEMBLE_SIZE\}/, model["ensemble_size_in_leaf"].to_s)
      param_str.gsub!(/\$\{COMP_RUNDIR\}/, "#{model[".run_path"]}/data")
      model_cmd.push(param_str)
    end
  end
  tempfile.puts Shellwords.escape(model_cmd.join(" "))
  parallel_settings = model["overridable_settings"]["parallel_settings"]
  tempfile.puts parallel_settings["num_total_procs"]
  tempfile.puts parallel_settings["num_threads"]
end
tempfile.close

status, stdout = ExternalScript.callScript("config/machines/#{case_configuration.machine}/submit", basic_env, true, tempfile.path)
logname = ScriptLogger.writeTextToFile("execute", stdout)

print "Finished. See \"#{logname}\" for execution log.\n"

File.unlink(tempfile.path)
