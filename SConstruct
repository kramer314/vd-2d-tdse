# Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
# See the LICENSE.txt file at the top-level directory of this distribution.

# Flag to specify custom build configuration
AddOption("--config-json", dest="config_json_fname", default="./config_json.default")


# Standard library imports
import os
import json


# Helper functions
def abspath(path):
    return os.path.abspath(path) + "/"


# Parse build configuration from JSON
config_json_fname = GetOption("config_json_fname")
config_json_f = open(config_json_fname, "r")
config_json = json.load(config_json_f)
config_json_f.close()


# Extract source and build directories
source_dir = abspath(config_json["src_path"])
build_dir = abspath(config_json["build_path"])


# Extract library names and paths
lib_names = []
lib_paths = []
lib_path_root = config_json["lib_path_root"]
libs_dict = config_json["libs"]

for lib in libs_dict:
  lib_names.append(lib["name"])

  lib_path = lib["path"]

  if (lib["use_root"]):
    lib_path = lib_path_root + lib_path

  lib_paths.append(abspath(lib_path))


# Extract compiler flags
flags = ""

debug = bool(config_json["debug"])

flags_dict = config_json["flags"]

for flag_cat, flag_params in flags_dict.items():
    flag_list = flag_params["list"]
    if (debug and bool(flag_params["debug"])) or (not debug and bool(flag_params["prod"])):
        flags += " ".join(flag for flag in flag_list) + " "


# Begin build
env = DefaultEnvironment(ENV = os.environ, TOOLS = ['default', "gfortran"])

env.Replace(F90FLAGS = flags)
env.Replace(LINKFLAGS = flags)
env.Replace(FORTRANMODDIRPREFIX = "-J ")
env.Replace(FORTRANMODDIR = build_dir)
env.Replace(F90PATH = [lib_paths, build_dir])

Export("env")
Export("lib_names")
Export("lib_paths")

SConscript(source_dir+"SConscript", variant_dir=build_dir, duplicate=1)


# For whatever reason, we can't use duplicate=0 and have *.mod files in the
# build directory. But, if we duplicate the source tree into the build
# directory SCons doesn't automatically clean the source files, so we have to
# manually define the entire build directory as a cleaning target.
Clean(".", build_dir)
