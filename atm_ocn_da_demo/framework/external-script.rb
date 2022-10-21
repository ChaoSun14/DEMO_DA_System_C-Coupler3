#!/usr/bin/env ruby

require "open3"
require "stringio"
require "shellwords"
require "tempfile"
require "pp"

module ExternalScript

  def self.callCommand(env, command, display_on_screen)
    stdout, status = Open3.popen2e(env, command) do
      |stdin, stdout, pipe_thread|

      semaphore = Mutex.new
      strout = StringIO.new

      stdin_thread = Thread.new do
        while !STDIN.eof? do
          begin
            str = STDIN.readpartial(1024)
            stdin.write(str)
            semaphore.synchronize {
              strout.write(str)
            }
          rescue Errno::EPIPE
            break
          end
        end
      end

      while !stdout.eof? do
        s = stdout.readpartial(1024)
        STDOUT.write(s) if display_on_screen
        semaphore.synchronize {
          strout.write(s)
        }
      end

      stdin_thread.terminate
      next strout.string, pipe_thread.value
    end
    return status, stdout
  end

  def self.callScript(script, env, display_on_screen, *args)
    tmpfile_list = Array.new
    cmdline = [ script ]
    args.each do |value|
      if value.is_a? Hash
        tmpfile = Tempfile.new(["testing", ".hash"])
        tmpfile_list.push(tmpfile)
        value.each do |hash_key, hash_value|
          tmpfile.puts("#{hash_key.upcase}=#{hash_value}")
        end
        tmpfile.close
        cmdline << Shellwords.escape(File.expand_path(tmpfile.path))
      elsif value.is_a? Array
        tmpfile = Tempfile.new(["testing", ".array"])
        tmpfile_list.push(tmpfile)
        value.each { |subvalue| tmpfile.puts(suvalue.to_s) }
        tmpfile.close
        cmdline << Shellwords.escape(File.expand_path(tmpfile.path))
      else
        cmdline << Shellwords.escape(value.to_s)
      end
    end
    return_tmpfile = Tempfile.new(["testing", ".return"])
    tmpfile_list.push(return_tmpfile)
    cmdline << Shellwords.escape(File.expand_path(return_tmpfile.path))

    status, stdout = callCommand(env, cmdline.join(" "), display_on_screen)

    return_text = return_tmpfile.read
    tmpfile_list.each{ |tmpfile| tmpfile.unlink }
    return status, stdout, return_text
  end

end
