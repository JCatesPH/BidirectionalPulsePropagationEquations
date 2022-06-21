#!/bin/bash
OUTPUT='test'
REMOVE='false'
FRAMERATE=30
PNGNAMES="frame%d.png"

print_help() {
    echo "Converts sequence of output%d.png files into mp4 using ffmpeg."
    echo
    echo "-r       : (default=false) sets whether to clean directory of png files"
    echo "-o OPT   : (default=test) sets the name of mp4 file with extension automatically added"
    echo "-f OPT   : (default=30) sets the frames per second of output video" 
    echo "-i REGEX : (default=frame) sets the regular expression for the frames \n\t[NOTE: numbers and file extension automatic]"
}

while getopts o:f:i:rh flag
do
    case "${flag}" in
        o) OUTPUT="${OPTARG}";;
        f) FRAMERATE="${OPTARG}";;
        i) PNGNAMES="${OPTARG}";;
        r) REMOVE='true';;
        h) print_help
            exit 1;;
    esac
done


# ffmpeg command to create mp4 of png figs
ffmpeg -r $FRAMERATE -f image2 -i "${PNGNAMES}%d.png" -vf scale=-2\:720 -vcodec libx264 -crf 25 -y -pix_fmt yuv420p "$OUTPUT".mp4

if [[ $REMOVE = 'true' ]]
then
    rm "${PNGNAMES}"*.png
fi
