#!/bin/bash
echo "Scipt location is ${0}" # full path to script including script name
echo -e " Initial wd is $(pwd), \n this is where you called the script from"
echo "Changing wd to...  ${0%/*}"
cd ${0%/*} # sets wd as location of the script, following symlinks
# cd "$(dirname "$0")" # same as above
# otherwise, wd is the pwd that was already set when script was called
echo "Now the wd is ... $(pwd)"
echo "with two arg passes, the wd would be ${0%/*}"