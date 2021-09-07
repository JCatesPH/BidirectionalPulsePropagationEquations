#!/bin/bash
OUTPUT='test'
REMOVE='false'
FRAMERATE=30

print_help() {
    echo "Converts sequence of output%d.png files into mp4 using ffmpeg."
    echo
    echo "-r     : (default=false) sets whether to clean directory of png files"
    echo "-o OPT : (default=test) sets the name of mp4 file with extension automatically added"
    echo "-f OPT : (default=30) sets the frames per second of output video" 
}

while getopts o:f:rh flag
do
    case "${flag}" in
        o) OUTPUT="${OPTARG}";;
        f) FRAMERATE="${OPTARG}";;
        r) REMOVE='true';;
        h) print_help
            exit 1;;
    esac
done


# ffmpeg command to create mp4 of png figs
ffmpeg -r $FRAMERATE -f image2 -i output%d.png -vf scale=-2:720 -vcodec libx264 -crf 25 -y -pix_fmt yuv420p "$OUTPUT".mp4

if [[ $REMOVE = 'true' ]]
then
    rm output*.png
fi
