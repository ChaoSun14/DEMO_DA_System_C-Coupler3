#!/usr/bin/env ruby

require "open3"
require "shellwords"

module GitOperator

  def self.isReallyClean(code_directory)
    ret, ret_msg = executeGitCommand(code_directory, "status", "--porcelain")
    if ret != 0
      return nil
    elsif ret_msg.length == 0
      return true
    else
      return false
    end
  end

  def self.getCurrentDigest(code_directory)
    ret, sha = executeGitCommand(code_directory, "rev-parse", "--verify", "HEAD")
    if ret == 0
      return sha.strip
    else
      return nil
    end
  end

  def self.getRefNameByDigest(code_directory, sha)
    ret, ref_name = executeGitCommand(code_directory, "name-rev", "--name-only", sha)
    if ret == 0
      return ref_name.strip
    else
      return nil
    end
  end

  def self.getDigestByTag(code_directory, tag)
    ret, sha = executeGitCommand(code_directory, "rev-parse", tag + "^{}")
    if ret == 0
      return sha.strip
    else
      return nil
    end
  end

  def self.getPatch(code_directory)
    ret, patch = executeGitCommand(code_directory, "diff")
    return nil if ret != 0
    ret, list = executeGitCommand(code_directory, "ls-files",
                                  "--others", "--exclude-standard")
    return nil if ret != 0
    list.split("\n").each() do |untracked_file|
      ret, this_patch = executeGitCommand(code_directory, "diff",
                                          "--", "/dev/null", untracked_file)
      patch << this_patch if ret == 0
    end
    patch
  end

  private
  def self.executeGitCommand(code_directory, *args)
    cmd = [ "git" ].concat(args)
    current_dir = File.expand_path(Dir.pwd)
    Dir.chdir(code_directory)
    begin
      ret_msg, status = Open3.capture2e(Shellwords.join(cmd))
      ret = status.exitstatus
    rescue Errno::ENOENT
      ret = -1
      ret_msg = "Git not found"
    end
    Dir.chdir(current_dir)
    return [ ret, ret_msg ]
  end

end
