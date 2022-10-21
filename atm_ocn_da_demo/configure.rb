#!/usr/bin/env ruby

framework_root_directory = File.absolute_path(File.dirname __FILE__)
$LOAD_PATH << "#{framework_root_directory}"

require "fileutils"
require "framework/xmlparser"
require "framework/external-script"
require "framework/environment"
require "framework/logger"
require "framework/git"
require "getoptlong"
require "tmpdir"
require "pp"

ScriptLogger.newLog("configure.log", newdir = true)

# Parse paramters
opts = GetoptLong.new(
  ['--reproducibility-archive', '-r', GetoptLong::NO_ARGUMENT],
  ['--reproducibility-check', '-t', GetoptLong::NO_ARGUMENT],
  ['--quick-digest', '-q', GetoptLong::NO_ARGUMENT]
)

flag_repro_archive = false
flag_repro_check = false
flag_quick_digest = false

opts.each do |opt, arg|
  case opt
  when '--reproducibility-archive'
    flag_repro_archive = true
  when '--reproducibility-check'
    flag_repro_check = true
    flag_repro_archive = false
  when '--quick-digest'
    flag_quick_digest = true
  end
end

# Load configuration
case_configuration = XmlParser::CaseConfiguration.new("config/case.xml")
models = case_configuration.getModelsFlatList

ScriptLogger.log.info("#{models.length} models loaded.")
print "#{models.length} models will be configured.\n"

# Construct directories
models.each() do |model|
  path = File.absolute_path("run/" + model["type"] + "/" + model["name"])
  path << "/#{model["ensemble_id"]}" if model["ensemble_id"] != nil
  model[".run_path"] = path
  FileUtils.mkdir_p(path + "/data")
  ScriptLogger.log.info("Created directory \"#{path}/data\".")
  print "Created directory \"#{path}/data\".\n"
end

# Construct temporary directory for reproducibility
if flag_repro_archive then
  repro_directory = Dir.mktmpdir()
  repro_data_list = "#{repro_directory}/.reproducibility.data"
  puts "[DEBUG] reproducibility directory: #{repro_directory}"
end

# Set some basic environment
basic_env = Hash.new
basic_env["CASE_NAME"] = case_configuration.case_name
basic_env["PATH"] = "#{$framework_root_directory}/buildtools/utils:#{ENV["PATH"]}"
basic_env["CASE_ROOT"] = File.absolute_path(".")
basic_env["QUICK_DIGEST"] = "1" if flag_quick_digest
basic_env["DIGEST_CHECKING"] = "1" if flag_repro_check
basic_env["REPRODUCIBILITY_DATA_LIST"] = repro_data_list if flag_repro_archive

# Call pre-configure script
if case_configuration.pre_configure != nil
  print "Peform pre-configuration \"#{case_configuration.pre_configure}\".."
  status, stdout, = ExternalScript.callScript(case_configuration.pre_configure, basic_env, false)
  logname = ScriptLogger.writeTextToFile("pre-configure", stdout)
  ScriptLogger.log.info("Pre-configure \"#{case_configuration.pre_configure}\", log: #{logname}")
  if status != 0 then
    print "error.\n\n"
    print "See \"#{logname}\" for more details.\n"
    exit
  end
  print "done.\n"
end

# Call configuration scripts of models
models.each() do |model|
  new_env, model_configuration = EnvironmentSetting.setEnvironment(basic_env, model)
  model["model_dir"] = model_configuration.model_dir
  print "Configure #{model["name"]}"
  print "(ensemble id ##{model["ensemble_id"]})" if model["ensemble_id"] != nil
  print ".. "
  status, stdout, = ExternalScript.callScript(model_configuration.configure_script, new_env, model_configuration.configure_interaction)
  logname = ScriptLogger.writeTextToFile("#{model["name"]}.config", stdout)
  ScriptLogger.log.info("Called configuration script of #{model["name"]}, log: #{logname}")
  if status != 0 then
    print "error.\n\n"
    print "See \"#{logname}\" for more details.\n"
    exit
  end
  print "done.\n"
end

# Call post-configure script
if case_configuration.post_configure != nil
  print "Peform post-configuration \"#{case_configuration.post_configure}\".."
  status, stdout, = ExternalScript.callScript(case_configuration.post_configure, basic_env, false)
  logname = ScriptLogger.writeTextToFile("post-configure", stdout)
  ScriptLogger.log.info("Post-configure \"#{case_configuration.post_configure}\", log: #{logname}")
  if status != 0 then
    print "error.\n\n"
    print "See \"#{logname}\" for more details.\n"
    exit
  end
  print "done."
end

if flag_repro_archive then
  code_info_dir = "#{repro_directory}/.reproducibility.code"
  FileUtils.mkdir_p(code_info_dir)
  model_exist = Hash.new
  models.each() do |model|
    model_name = "#{model["name"]}"
    if model_exist.has_key? model_name then
      next
    else
      model_exist[model_name] = true
    end
    model_dir = model["model_dir"]
    digest = GitOperator.getCurrentDigest(model_dir)
    is_clean = GitOperator.isReallyClean(model_dir)
    if ! is_clean then
      patch = GitOperator.getPatch(model_dir)
      File.open("#{code_info_dir}/#{model_name}.patch", "w") do |f|
        f.write(patch)
      end
    end
    File.open("#{code_info_dir}/version.txt", "a") do |f|
      f.write("#{model_name} : #{digest}")
      if is_clean then
        f.write("\n")
      else
        f.write(", #{model_name}.patch\n")
      end
    end
  end
  FileUtils.cp_r("config", repro_directory)
  if Dir.exist?("CCPL_dir") then
    FileUtils.mkdir_p("#{repro_directory}/CCPL_dir")
    FileUtils.mkdir_p("#{repro_directory}/CCPL_dir/run")
    FileUtils.cp_r("CCPL_dir/config", "#{repro_directory}/CCPL_dir/")
    FileUtils.cp_r("CCPL_dir/run/data", "#{repro_directory}/CCPL_dir/run/")
  end

  FileUtils.mkdir_p("reproducibility")
  archive_file = File.absolute_path("reproducibility/#{ScriptLogger.nowDateTimeStr}.tar")
  Dir.chdir(repro_directory) do
    system("tar cf #{archive_file} * .reproducibility.*")
    FileUtils.rm_rf repro_directory
  end

end

print "Configuration finished.\n"
