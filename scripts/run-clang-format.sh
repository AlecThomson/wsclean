
#!/bin/bash
#
# Author: Jakob Maljaars
# Email: jakob.maljaars_@_stcorp.nl

programname=$0

set -e

# Function to print usage of this script
usage () {
    s=$(printf "%-50s" "=")
    echo "${s// /=}"
    echo "Script to format C++ header and source files. Expects a .clang-format file in current directory or parent"
    echo ""
    echo "Usage: $programname [-i filename/pattern] [-s skip directory]"
    echo "  -i      include file(pattern) in diff"
    echo "  -s      skip directory in diff. Must be specified as a path, e.g. ./external"
    echo "  -h      display help"
    exit 1
}

# Function to join array to a string
join_by () { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

# Fill the include_files and exclude_directories
include_files=()
exclude_directories=()
while getopts ":hi:s:" opt
do
    case $opt in
        h) 
        usage
        exit;;
        i) 
        include_files=("${include_files[@]}" $OPTARG)
        ;;
        s) 
        exclude_directories=("${exclude_directories[@]}" $OPTARG)
        ;;
        \? ) 
        echo "Execute $programname -h to get usage of this script."
        exit;;
    esac
done

# Provide a default set of file extensions if -i flag not set
if [ ${#include_files[@]} -eq 0 ]; then
    file_ext=(*.c *.cc *.cxx *.cpp *.c++ *.hh *.hxx *.hpp *.h)
else
    file_ext=("${include_files[@]}")
fi
# Join into variable
INC_FILE=$(join_by " -o -name " "${file_ext[@]}")

if [ ${#exclude_directories[@]} -eq 0 ]; then
    EX_DIR=""
else
    EX_DIR=$(join_by " -path " "${exclude_directories[@]}")
fi

find . -path $EX_DIR -prune -o -type f \( -name $INC_FILE \) -exec clang-format -i -style=file \{\} +