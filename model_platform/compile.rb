#!/usr/bin/env ruby

framework_root_directory = File.absolute_path(File.dirname __FILE__)
$LOAD_PATH << "#{framework_root_directory}"
$framework_root_directory = framework_root_directory

require "fileutils"
require "framework/xmlparser"
require "framework/external-script"
require "framework/environment"
require "framework/logger"
require "open3"
require "digest"
require "shellwords"
require "tempfile"
require "pp"

def createDstDirs(dst_path)
  FileUtils.mkdir_p("#{dst_path}/src")
  FileUtils.mkdir_p("#{dst_path}/obj")
  FileUtils.mkdir_p("#{dst_path}/exe")
  ScriptLogger.log.info("Created directory \"#{dst_path}/src\".")
  ScriptLogger.log.info("Created directory \"#{dst_path}/obj\".")
  ScriptLogger.log.info("Created directory \"#{dst_path}/exe\".")
  print "Created directory \"#{dst_path}/src\"\n"
  print "Created directory \"#{dst_path}/obj\"\n"
  print "Created directory \"#{dst_path}/exe\"\n"
end

def defaultCompile(src_path, dst_path, exefile, env, pre_build, form_src, is_library = false)
  dst_path_src = dst_path + "/src"
  dst_path_obj = dst_path + "/obj"
  status, stdout, return_text = ExternalScript.callScript(form_src, env, false, src_path)
  srclist = return_text.split("\n")
  srclist.each do |srcpath|
    raise StandardError.new "#{srcpath} does not exist." if ! File.exist?(srcpath)
    Dir.chdir(srcpath) do
      Dir.glob("*.{h,c,cxx,f,F,f90,F90}") { |file|
        FileUtils.cp_r(file, dst_path_src, preserve: true)
      }
    end
    status, stdout, return_text = ExternalScript.callScript(pre_build, env, false) if File.exist?(pre_build)
    Dir.chdir(dst_path_src) do
      system("grep  \"^[\\t ]*module \" -Prni * | grep -vni \"procedure\" | sed \"s/:/    /g\"  |  awk '{print $2 \" \" $5}' > #{dst_path_obj}/module_mapping.info")
    end
    system("cc #{$framework_root_directory}/buildtools/makdep.c -o #{dst_path_obj}/makdep")
  end
  Dir.chdir(dst_path_obj) do
    if is_library then
      cmdline = "gmake -j16 -f #{$framework_root_directory}/buildtools/Makefile.libs"
    else
      cmdline = "gmake -j16 -f #{$framework_root_directory}/buildtools/Makefile.models"
    end
    new_env = env.clone
    File.write("Srcfiles", "#{dst_path_src}\n")
    #puts(File.new(new_env["MACFILE"]).read)
    new_env["VPFILE"] = File.absolute_path("Srcfiles")
    new_env["EXEC"] = exefile
    status, stdout = ExternalScript.callCommand(new_env, cmdline, false)
  end
  return status, stdout
end

def generateCompileCfg(common_compiler_cfg, component_compiler_cfg,
                       config, config_dir, libraries = nil)
  cfg_list = Array.new
  file = Tempfile.new("compiler_cfg_")
  file.write(File.read(common_compiler_cfg))
  cfg_list.push(file)
  file.close

  if component_compiler_cfg != nil \
      and File.exist?(component_compiler_cfg) then
    text = File.read(component_compiler_cfg)
    text_new = text.gsub(/#PARENT_COMPILER_CFG#/, cfg_list.last.path)
    file = Tempfile.new("compiler_cfg_")
    cfg_list.push(file)
    file.write(text_new)
    file.close
  end

  if config["overridable_settings"]["compiling_settings"] != nil \
    and ! config["overridable_settings"]["compiling_settings"].empty? then
    compiling_settings =
      config["overridable_settings"]["compiling_settings"]
    if compiling_settings.has_key?("file") then
      filename = compiling_settings["file"]
      filename = "#{config_dir}/#{filename}" if filename[0] != "/"
      if File.exist?(filename) then
        text = File.read(filename)
        text_new = text.gsub(/#PARENT_COMPILER_CFG#/, cfg_list.last.path)
        tmpfile = Tempfile.new("compiler_cfg_")
        cfg_list.push(tmpfile)
        tmpfile.write(text_new)
        tmpfile.close
      else
        raise "Compiler configuration file \"#{compiling_settings["file"]}\" was not found."
      end
    end
  end

  if libraries != nil then
    ldflags_string = ""
    include_string = ""
    libraries.each do |library|
      ldflags_string += " -l#{library["name"]} -L#{library["path"]}"
      include_string += " -I#{library["include"]}"
    end
    file = Tempfile.new("compiler_cfg_")
    file.write("include #{cfg_list.last.path}\n")
    file.write("SLIBS += #{ldflags_string}\n")
    file.write("INCLDIR += #{include_string}\n")
    cfg_list.push(file)
    file.close
  end
  return cfg_list
end

def buildLibrary(library, env, common_compiler_cfg, library_hash)
  library_path = File.absolute_path(
    "run/libs/" + library["name"] + "/" +
    Digest::SHA1.hexdigest(Marshal.dump(library)))
  FileUtils.mkdir_p(library_path)
  library_config_path = File.absolute_path(
    "config/libs/" + library["name"])
  library_config = XmlParser::LibraryConfiguration.new(library_config_path + "/lib.xml", library["type"], library["name"])
  cfg_list = generateCompileCfg(common_compiler_cfg, "#{library_config_path}/compiler.cfg", library, library_config)
  env["MACFILE"] = cfg_list.last.path
  createDstDirs(library_path)
  print "Compile library \"#{library["name"]}\" "
  if library_config.use_default_compile_script then
    print "with default compiling script .."
    status, stdout =
      defaultCompile(library_config.library_dir, library_path,
                     library_path + "/exe/lib#{library["name"]}.a", env,
                     library_config_path + "/pre_build.sh",
                     library_config_path + "/form_src.sh", true)
    logname = ScriptLogger.writeTextToFile("#{library["name"]}.#{library_hash}.compile", stdout)
    ScriptLogger.log.info("Compiled library \"#{library["name"]}\" with default compiling script, log: \"#{logname}\"")
  else
    print "with custom compiling script \"#{library_config.compile_script}\".. "
    status, stdout, = ExternalScript.callScript(
      library_config.compile_script, new_env, false)
    logname = ScriptLogger.writeTextToFile("#{library["name"]}.#{library_hash}.compile", stdout)
    ScriptLogger.log.info("Compiled library \"#{library["name"]}\" with custom compiling script \"#{library_config.compile_script}\", log: \"#{logname}\"")
  end
  cfg_list.each do |cfg_file|
    cfg_file.unlink
  end
  if status != 0 then
    print "error.\n\n"
    printf "See \"#{logname}\" for more details.\n"
    exit
  end
  print "done.\n"
  return library["name"], "#{library_path}/exe", "#{library_path}/obj"
end

def buildModel(model, env, common_compiler_cfg, libraries)
  model_path = model[".run_path"]
  model_config_path = File.absolute_path(
    "config/models/#{model["type"]}/#{model["name"]}")
  model_config = XmlParser::ModelConfiguration.new(
    "#{model_config_path}/model.xml", model["type"], model["name"])
  cfg_list = generateCompileCfg(common_compiler_cfg, "#{model_config_path}/compiler.cfg", model, model_config_path, libraries)
  env["MACFILE"] = cfg_list.last.path
  createDstDirs(model_path)
  print "  Compile model \"#{model["name"]}\" "
  logfile_prefix = model["name"]
  logfile_prefix = "#{logfile_prefix}.#{model["ensemble_id"]}" if model["ensemble_id"] != nil
  if model_config.use_default_compile_script then
    print "with default compiling script "
    status, stdout =
      defaultCompile(model_config.model_dir, model_path,
                     "#{model_path}/exe/#{model["name"]}", env,
                     "#{model_config_path}/pre_build.sh",
                     "#{model_config_path}/form_src.sh", false)
    logname = ScriptLogger.writeTextToFile("#{logfile_prefix}.compile", stdout)
    ScriptLogger.log.info("Compiled model \"#{logfile_prefix}\" with default compiling script, log: \"#{logname}\"")
  else
    print "with custom compiling script \"#{model_config.compile_script}\" "
    env["EXEFILE"] = "#{model_path}/exe/#{model["name"]}"
    status, stdout,  = ExternalScript.callScript(model_config.compile_script, env, true)
    logname = ScriptLogger.writeTextToFile("#{logfile_prefix}.compile", stdout)
    ScriptLogger.log.info("Compiled model \"#{logfile_prefix}\" with custom compiling script \"#{model_config.compile_script}\", log: \"#{logname}\"")
  end
  cfg_list.each do |cfg_file|
    cfg_file.unlink
  end
  if status != 0 then
    print "error.\n\n"
    printf "See \"#{logname}\" for more details.\n"
    exit
  end
  print "done.\n"
end

def copyModelExec(model, execfile)
  dst_execfile = "#{model[".run_path"]}/exe/#{model[".name"]}"
  createDstDirs(model[".run_path"])
  FileUtils.cp(execfile, dst_execfile)
end

ScriptLogger.newLog("compile.log", newdir = false, append_datetime = true)

# Parse paramters

# Load configuration
case_configuration = XmlParser::CaseConfiguration.new("config/case.xml")
models = case_configuration.getModelsFlatList

models.each() do |model|
  path = File.absolute_path("run/" + model["type"] + "/" + model["name"])
  path << "/#{model["ensemble_id"]}" if model["ensemble_id"] != nil
  model[".run_path"] = path
end

ScriptLogger.log.info("#{models.length} models loaded.")
print "#{models.length} models will be compiled.\n"

# Set some basic environment
basic_env = Hash.new
basic_env["CASE_NAME"] = case_configuration.case_name
basic_env["PATH"] = "#{$framework_root_directory}/buildtools/utils:#{ENV["PATH"]}"
basic_env["CASE_ROOT"] = File.absolute_path(".")

common_compiler_cfg = File.absolute_path("config/machines/#{case_configuration.machine}/common_compiler.cfg")

models_has_been_built = Hash.new
libraries_has_been_built = Hash.new
models.each() do |model|
  new_env, model_configuration = EnvironmentSetting.setEnvironment(basic_env, model)
  libraries_settings = model["overridable_settings"]["libraries"] ||
    Array.new
  libraries_settings_for_makefile = Array.new
  libraries_settings.each() do |library|
    library_hash = Digest::SHA1.hexdigest(Marshal.dump(library))
    if ! libraries_has_been_built[library_hash] then
      library_name, library_path, include_path =
        buildLibrary(library, new_env, common_compiler_cfg, library_hash)
      libraries_has_been_built[library_hash] = {"name" => library_name, "path" => library_path, "include" => include_path}
    end
    libraries_settings_for_makefile.push(
      libraries_has_been_built[library_hash])
  end
  if model["rebuild"] == "NO" then
    model_tmp = Marshal.load(Marshal.dump(model["default_settings"]))
  else
    model_tmp = Marshal.load(Marshal.dump(model))
  end
  model_tmp.delete("ensemble_id")
  model_tmp.delete("ensemble_size_in_leaf")
  model_tmp.delete(".run_path")
  model_tmp.delete("rebuild")
  model_key = Digest::SHA1.hexdigest(Marshal.dump(model_tmp))

  print "Compile #{model["name"]}"
  print "(ensemble ##{model["ensemble_id"]}) " if model["ensemble_id"] != nil
  if model["rebuild"] == "YES" then
    print ": forced rebuilding)..\n"
    buildModel(model, new_env, common_compiler_cfg, libraries_settings_for_makefile)
  elsif models_has_been_built.has_key?(model_key) then
    print ": dectected built executable file.\n"
    copyModelExec(model, models_has_been_built[model_key])
  elsif model["rebuild"] == "NO" then
    print ": force not rebuilding, but default executable does not exist. Will build it...\n"
    model_tmp["ensemble_id"] = model["ensemble_id"]
    model_tmp[".run_path"] = model[".run_path"]
    buildModel(model_tmp, new_env, common_compiler_cfg, libraries_settings_for_makefile)
    models_has_been_built[model_key] = "#{model_tmp[".run_path"]}/exe/#{model_tmp["name"]}"
  else
    print ": model has not been built, will build it...\n"
    buildModel(model, new_env, common_compiler_cfg, libraries_settings_for_makefile)
    models_has_been_built[model_key] = "#{model[".run_path"]}/exe/#{model["name"]}"
  end
end

print "Compilation finished.\n"
