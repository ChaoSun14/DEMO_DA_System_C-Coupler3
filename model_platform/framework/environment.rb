#!/usr/bin/env ruby

require "framework/xmlparser"

module EnvironmentSetting

  def self.setEnvironment(old_env, model_config)
    new_env = old_env.clone

    new_env["RUN_PATH"] = model_config[".run_path"]
    new_env["MODEL_TYPE"] = model_config["type"]
    new_env["MODEL_NAME"] = model_config["name"]
    new_env["DATA_DIR"] = new_env["RUN_PATH"] + "/data"
    new_env["ENSEMBLE_ID"] = model_config["ensemble_id"].to_s

    running_settings = model_config["overridable_settings"]["running_settings"]
    new_env["RUN_TYPE"] = running_settings["run.type"]
    running_settings["run"].each do |key, value|
      new_env["RUNNING_SETTINGS_RUN_" + key.upcase] = value
    end
    running_settings["write_restart"].each do |key, value|
      new_env["RUNNING_SETTINGS_WRITE_RESTART_" + key.upcase] = value
    end
    parallel_settings = model_config["overridable_settings"]["parallel_settings"]
    new_env["PARALLEL_SETTINGS_NUM_TOTAL_PROCS"] = parallel_settings["num_total_procs"].to_s
    new_env["PARALLEL_SETTINGS_NUM_THREADS"] = parallel_settings["num_threads"].to_s
    parallel_settings["custom_settings"].each do |key, value|
      new_env["PARALLEL_SETTINGS_" + key.upcase] = value
    end
    custom_settings = model_config["overridable_settings"]["custom_settings"]
    custom_settings.each do |key, value|
      new_env["CUSTOM_SETTINGS_" + key.upcase] = value
    end
    configure_path = "config/models/" + model_config["type"] + "/" + model_config["name"]
    model_configuration = XmlParser::ModelConfiguration.new(configure_path + "/model.xml", model_config["type"], model_config["name"])
    new_env["MODEL_DATA_DIR"] = model_configuration.data_dir
    new_env["MODEL_DIR"] = model_configuration.model_dir
    return new_env, model_configuration
  end

end
