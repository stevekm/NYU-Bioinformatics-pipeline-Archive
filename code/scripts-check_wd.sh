#!/bin/bash
echo " $(basename ${0%/*/*} ) " # <<<---- THIS is the one I want # return just the project dir, two levels up from the script
echo " ${0%/*/*} " # file path that script was called from, two dirs up
cd ${0%/*/*}
echo $(pwd)
