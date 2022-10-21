#!/usr/bin/env ruby

require "rexml/document"
require "set"
require "pp"

module XmlParser

  class BaseXmlConfiguration

    def initialize(xml_file)
      file = File.new(xml_file)
      @xml_doc = REXML::Document.new(file)
    end

    def getXmlText(xml_object, path, default_value_when_not_found = nil)
      ret_value = ""
      xml_object.elements.each(path) do |e|
        ret_value = e.text
      end.empty? and \
      begin
        ret_value = default_value_when_not_found
      end
      return ret_value
    end

    def parseChildrenToHash(xml_object, original_hash = nil)
      return original_hash if xml_object == nil
      hash_object = original_hash
      hash_object = Hash.new if hash_object == nil
      xml_object.elements.each() do |child|
        hash_object[child.name] = child.text
      end
      return hash_object
    end

    protected :getXmlText, :parseChildrenToHash

  end

  class CaseConfiguration < BaseXmlConfiguration

    def initialize(xml_file)
      super(xml_file)
      @case_name = getXmlText(@xml_doc, "case/name", "")
      @machine = getXmlText(@xml_doc, "case/machine", "generic_linux")
      @ccpl_default_namelist = getXmlText(@xml_doc, "case/ccpl_default_namelist", "yes")
      @pre_configure = getXmlText(@xml_doc, "case/pre_configure")
      @post_configure = getXmlText(@xml_doc, "case/post_configure")
      @machine_options = 
        parseChildrenToHash(@xml_doc.root.elements["machine_options"])
      overridable_settings = Hash.new
      overridable_settings["custom_settings"] =
        parseChildrenToHash(@xml_doc.root.elements["custom_settings"])
      overridable_settings["running_settings"] =
        parseRunningSettings(@xml_doc.root.elements["running_settings"])
      overridable_settings["compiling_settings"] =
        parseCompilingSettings(@xml_doc.root.elements["compiling_settings"])
      overridable_settings["submit_settings"] =
        parseSubmitSettings(@xml_doc.root.elements["submit_settings"])
      overridable_settings["parallel_settings"] =
        parseParallelSettings(@xml_doc.root.elements["parallel_settings"])
      @libraries = parseLibraries(@xml_doc.root.elements["libraries"],
                                  overridable_settings)
      overridable_settings["libraries"] = @libraries
      @models = parseModels(@xml_doc.root.elements["models"],
                            overridable_settings)
    end

    def deep_copy(o)
      Marshal.load(Marshal.dump(o))
    end

    def parseRunningSettings(xml_object, original_settings = nil)
      return original_settings if xml_object == nil
      running_settings = original_settings
      running_settings = Hash.new if running_settings == nil
      run_node = xml_object.elements["run"]
      running_settings["run.type"] = run_node.attributes["type"]
      running_settings["run"] =
        parseChildrenToHash(run_node, running_settings["run"])
      running_settings["write_restart"] =
        parseChildrenToHash(xml_object.elements["write_restart"],
                            running_settings["write_restart"])
      return running_settings
    end

    def parseCompilingSettings(xml_object, original_settings = nil)
      return original_settings if xml_object == nil
      compiling_settings = original_settings
      compiling_settings = Hash.new if compiling_settings == nil
      if xml_object.elements["file"] != nil then
        compiling_settings["file"] =
          xml_object.elements["file"].attributes["filename"]
      end
      xml_object.elements.each("flag") do |child|
        compiling_settings["flag"] = Array.new if !compiling_settings.has_key?("flag")
        compiling_settings["flag"].push(child.attributes.to_a)
      end
      return compiling_settings
    end

    def parseSubmitSettings(xml_object, original_settings = nil)
      return original_settings if xml_object == nil
      submit_settings = Array.new
      xml_object.elements.each("paramter") do |child|
        submit_settings.push(child.text)
      end
      return submit_settings
    end

    def parseParallelSettings(xml_object, original_settings = nil)
      return original_settings if xml_object == nil
      parallel_settings = Hash.new
      parallel_settings["num_total_procs"] =
        (xml_object.attributes["num_total_procs"] || "1").to_i
      parallel_settings["num_threads"] =
        (xml_object.attributes["num_threads"] || "1").to_i
      parallel_settings["custom_settings"] =
        parseChildrenToHash(xml_object)
      return parallel_settings
    end

    def parseLibraries(xml_object, overridable_settings,
                       original_libraries = nil)
      return original_libraries if xml_object == nil
      is_append = (xml_object.attributes["append"] != "no")
      libraries = original_libraries
      libraries = Array.new if libraries == nil || ! is_append
      xml_object.elements.each("library") do |child|
        this_library = Hash.new
        this_library["type"] = child.attributes["type"]
        this_library["name"] = child.attributes["name"]
        settings = Hash.new
        settings["compiling_settings"] =
          parseCompilingSettings(child.elements["compiling_settings"],
                                 overridable_settings["compiling_settings"]
                                )
        settings["custom_settings"] =
          parseChildrenToHash(child.elements["custom_settings"],
                              overridable_settings["custom_settings"]
                             )
        this_library["overridable_settings"] = deep_copy(settings)
        libraries.push(this_library)
      end
      return libraries
    end

    def parseAllOverridableSettings(xml_object, original_settings = nil)
      return original_settings if xml_object == nil
      settings = deep_copy(original_settings) || Hash.new
      settings["custom_settings"] = parseChildrenToHash(
        xml_object.elements["custom_settings"],
        settings["custom_settings"]
      )
      settings["running_settings"] = parseRunningSettings(
        xml_object.elements["running_settings"],
        settings["running_settings"]
      )
      settings["compiling_settings"] = parseCompilingSettings(
        xml_object.elements["compiling_settings"],
        settings["compiling_settings"]
      )
      settings["submit_settings"] = parseSubmitSettings(
        xml_object.elements["submit_settings"],
        settings["submit_settings"]
      )
      settings["parallel_settings"] = parseParallelSettings(
        xml_object.elements["parallel_settings"],
        settings["parallel_settings"]
      )
      settings["libraries"] = parseLibraries(
        xml_object.elements["libraries"],
        settings, settings["libraries"]
      )
      return settings
    end

    def parseId(id_string)
      id_set = Set.new
      id_string.split(",").each() do |range|
        range_array = range.split("-")
        if range_array.count == 1 then
          id_set.add(range_array[0].to_i)
        else
          id_set.merge((range_array[0].to_i)..(range_array[1].to_i))
        end
      end
      return id_set
    end

    def parseModel(xml_object, overridable_settings, in_ensemble = false,
                  in_coupled = false, member_template = nil)
      return nil if xml_object == nil
      this_model = member_template || Hash.new
      this_model["type"] = xml_object.attributes["type"] ||
        this_model["type"]
      this_model["name"] = xml_object.attributes["name"] ||
        this_model["name"]
      this_model["ensemble_size"] =
        (xml_object.attributes["ensemble_size"] || "0").to_i
      this_model["code_version"] =
        xml_object.attributes["code_version"] ||
        this_model["code_version"] || "-"
      is_coupled = this_model["type"] == "coupled"
      model_settings_node = xml_object
      model_settings_node = xml_object.elements["default"] if this_model["ensemble_size"] > 1
      settings = parseAllOverridableSettings(model_settings_node,
                                             overridable_settings)
      this_model["overridable_settings"] = settings
      this_model["..is_leaf"] = true
      if is_coupled then
        this_model["models"] = parseModels(
          model_settings_node.elements["models"],
          settings, false, true)
        this_model["..is_leaf"] = false
      end
      if this_model["ensemble_size"] > 1 then
        this_model["..is_leaf"] = false
        ensemble_members = Array.new
        for i in 1..this_model["ensemble_size"] do
          member = deep_copy(this_model)
          member["ensemble_size_in_leaf"] = member["ensemble_size"]
          member.delete("ensemble_size")
          member["ensemble_id"] = i
          member["..is_leaf"] = true
          ensemble_members.push(member)
        end
        specific_members_node =
          xml_object.elements["specific_ensenble_members"]
        if specific_members_node != nil then
          specific_members_node.elements.each("member") do |child|
            next if child.attributes["status"] == "off"
            id_set = parseId(child.attributes["id"])
            for id in id_set do
              if is_coupled then
                ensemble_members[id - 1] =
                  parseModels(child.elements["models"], settings, true,
                              true, ensemble_members[id - 1])
              else
                ensemble_members[id - 1] =
                  parseModel(child, settings, true, false,
                            ensemble_members[id - 1])
              end
              ensemble_members[id - 1]["rebuild"] = (child.attributes["rebuild"] || "auto").upcase
              if ensemble_members[id - 1]["rebuild"] == "NO" then
                ensemble_members[id - 1]["default_settings"] = deep_copy(this_model)
              end
            end
          end
        end
        this_model["ensemble_members"] = ensemble_members
      end
      return this_model
    end

    def parseModels(xml_object, overridable_settings, in_ensemble = false,
                    in_coupled = false, member_template = nil)
      return nil if xml_object == nil
      if in_ensemble && xml_object.attributes["override"] == "none" then
        return member_template
      elsif in_ensemble && xml_object.attributes["override"] == "part" then
        models = member_template
        xml_object.elements.each("model") do |child|
          type = child.attributes["type"]
          name = child.attributes["name"]
          is_new_model = true
          models["models"].each_index do |i|
            model = models["models"][i]
            if model["type"] == type && model["name"] == name then
              if child.attributes["delete"] == "yes" then
                models["models"].delete(model)
              else
                models["models"][i] = 
                  parseModel(child, overridable_settings, true, true,
                                   model)
              end
              is_new_model = false
              break
            end
          end
          models["models"].push(child, overridable_settings,
                                true, true, nil) if is_new_model
        end
        return models
      else
        models = Array.new
        xml_object.elements.each("model") do |child|
          this_model = parseModel(child, overridable_settings, in_ensemble,
                                  in_coupled, member_template)
          models.push(this_model)
        end
        return models
      end
    end

    def getLeafModels(model, parent = nil)
      model_list = Array.new
      if (model["ensemble_size"] || 0) > 1 then
        model["ensemble_members"].each do |member|
          model_list.concat(getLeafModels(member))
        end
      elsif model["type"] == "coupled" then
        model["models"].each do |submodel|
          model_list.concat(getLeafModels(submodel, model))
        end
      else
        if model["..is_leaf"] then
          model_instance = deep_copy(model)
          if parent != nil then
            model_instance["ensemble_id"] = parent["ensemble_id"] if parent.has_key?("ensemble_id")
            model_instance["ensemble_size_in_leaf"] = parent["ensemble_size_in_leaf"] if parent.has_key?("ensemble_size_in_leaf")
            model_instance["rebuild"] = parent["rebuild"] if parent.has_key?("rebuild")
            model_instance["default_settings"] = parent["default_settings"] if parent.has_key?("default_settings")
            if parent.has_key?("default_settings") then
              parent["default_settings"]["models"].each do |submodel|
                if submodel["name"] == model["name"] &&
                    submodel["type"] == model["type"] then
                  model_instance["default_settings"] = submodel
                  break
                end
              end
            end
          end
          model_list.push(model_instance)
        end
      end
      model_list
    end

    def getModelsFlatList
      flat_list = Array.new
      @models.each do |model|
        flat_list.concat(getLeafModels(model))
      end
      flat_list
    end

    attr_reader :case_name, :ccpl_default_namelist, \
      :pre_configure, :post_configure, :libraries, \
      :models, :machine, :machine_options

    protected :parseRunningSettings, :parseCompilingSettings, \
      :parseSubmitSettings, :parseChildrenToHash, \
      :parseParallelSettings, :parseAllOverridableSettings, \
      :parseLibraries, :parseModel, :parseModels, \
      :getLeafModels
  end

  class ModelConfiguration < BaseXmlConfiguration
    def initialize(xml_file, model_type, model_name)
      super(xml_file)

      @model_dir = getXmlText(@xml_doc, "model/model_dir", nil)
      if @model_dir == nil then
        @model_dir = ENV["MODEL_ROOT"] + "/" + model_type + "/" + model_name
      elsif @model_dir[0] != "/" then
        @model_dir = ENV["MODEL_ROOT"] + "/" + @model_dir
      end
      @model_dir = File.absolute_path(@model_dir)

      @data_dir = getXmlText(@xml_doc, "model/data_dir", nil)
      if @data_dir == nil then
        @data_dir = ENV["DATA_ROOT"] + "/" + model_type + "/" + model_name
      elsif @data_dir[0] != "/" then
        @data_dir = ENV["DATA_ROOT"] + "/" + @data_dir
      end
      @data_dir = File.absolute_path(@data_dir)

      @configure_script = getXmlText(@xml_doc, "model/configure_script", nil)
      if @configure_script == nil then
        @configure_script = "config/models/" + model_type + "/" +
          model_name + "/config.sh"
      elsif @configure_script[0] != "/" then
        @configure_script = @model_dir + "/" + @configure_script
      end
      @configure_script = File.absolute_path(@configure_script)

      @compile_script = getXmlText(@xml_doc, "model/compile_script", nil)
      @use_default_compile_script = false
      if @compile_script == nil then
        @use_default_compile_script = true
      elsif @compile_script[0] != "/" then
        @compile_script = "config/models/#{model_type}/#{model_name}/#{@compile_script}"
      end
      @compile_script = File.absolute_path(@compile_script) if @compile_script != nil

      @configure_interaction = (getXmlText(@xml_doc, "model/configure_interaction", "no") == "yes")
    end

    attr_reader :model_dir, :data_dir, \
      :configure_script, :compile_script, \
      :use_default_compile_script, :configure_interaction
  end

  class LibraryConfiguration < BaseXmlConfiguration
    def initialize(xml_file, library_type, library_name)
      super(xml_file)

      @library_dir = getXmlText(@xml_doc, "library/library_dir", nil)
      if @library_dir == nil then
        @library_dir = ENV["MODEL_ROOT"] + "/libs/" + library_name
      elsif @library_dir[0] != "/" then
        @library_dir = ENV["MODEL_ROOT"] + "/" + @library_dir
      end
      @library_dir = File.absolute_path(@library_dir)

      @compile_script = getXmlText(@xml_doc, "library/compile_script", nil)
      @use_default_compile_script = false
      if @compile_script == nil then
        @use_default_compile_script = true
      elsif @compile_script[0] != "/" then
        @compile_script = "config/libs/#{library_name}/#{@compile_script}"
      end
      @compile_script = File.absolute_path(@compile_script) if @compile_script != nil

    end

    attr_reader :library_dir, :compile_script, :use_default_compile_script
  end

end
