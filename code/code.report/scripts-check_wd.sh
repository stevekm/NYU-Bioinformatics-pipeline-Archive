#!/bin/bash
echo "initial wd is :"
echo "$0"
echo " $(basename ${0%/*/*} ) " # return just the project dir, two levels up from the script
echo " ${0%/*/*} " # file path that script was called from, two dirs up
cd ${0%/*/*}

echo "The pwd is: "
echo $(pwd)
