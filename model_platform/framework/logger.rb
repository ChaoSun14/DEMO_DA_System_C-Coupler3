#!/usr/bin/env ruby

require "fileutils"
require "time"
require "logger"

module ScriptLogger

  @logger = nil
  @logdir = nil
  @now_date_time = nil

  def self.log
    if @logger.nil?
      @logger = Logger.new STDOUT
      @logger.level = Logger::DEBUG
      @logger.datetime_format = '%Y-%m-%d %H:%M:%S '
    end
    @logger
  end

  def self.newLog(filename, newdir = false, append_datetime = false)
    @logger.close if ! @logger.nil?
    if newdir then
      File.unlink(".last_configure_time") if File.exist?(".last_configure_time")
      getLogDirectory(true)
    end
    filename = "#{filename}.#{nowDateTimeStr}" if append_datetime
    @logger = Logger.new "#{getLogDirectory}/#{filename}"
    @logger.level = Logger::DEBUG
    @logger.datetime_format = '%Y-%m-%d %H:%M:%S '
    @logger
  end

  def self.getLogDirectory(renew = false)
    return @logdir if ! renew and ! @logdir.nil?
    if File.exist? ".last_configure_time" then
      time_str = File.open(".last_configure_time").readline.strip
    else
      time_str = nowDateTimeStr
      timefile = File.open(".last_configure_time", mode = "w")
      timefile.puts(time_str)
      timefile.close
    end
    @logdir = File.absolute_path("logs/#{time_str}")
    FileUtils.mkdir_p(@logdir)
    @logdir
  end

  def self.writeTextToFile(filename_prefix, text)
    filename = "#{getLogDirectory}/#{filename_prefix}.#{nowDateTimeStr}"
    fw = File.open(filename, "w")
    fw.write(text)
    fw.close
    filename
  end

  def self.nowDateTimeStr
    @now_date_time = Time.now.to_datetime.strftime("%Y%m%d-%H%M%S%z") if @now_date_time == nil
    @now_date_time
  end

end
