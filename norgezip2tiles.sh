#!/usr/bin/env bash

# Process a Geonorge N50 raster
# This script takes a geotiff from the goenorge database and converts
# it to web tiles, uploading to an s3 bucket.

tile_script=$1
zip_file=$2
PROJINFO="epsg:3857"
XWIDTH=4.77732
YWIDTH=4.77732
RESAMPLING="cubic"
SERVER="http://s3-us-west-1.amazonaws.com/caltopodev/norway-test/"
BUCKET="caltopodev"

tif_filename=$(unzip -l $zip_file | grep -oE "[^ ]+\.tif")
dir_name=$(basename $tif_filename .tif)_dir
tif_file="$dir_name/$tif_filename"
extension="${tif_filename##*.}"
name="${tif_filename%.*}"
warpedfile="${dir_name}/${name}_warped.${extension}"

mkdir $dir_name
unzip -d $dir_name $zip_file $tif_filename
gdalwarp -t_srs $PROJINFO -dstalpha -tr $XWIDTH $YWIDTH -r $RESAMPLING $tif_file $warpedfile
$tile_script -z 5-15 -r bilinear -b $BUCKET -f $SERVER $warpedfile $dir_name

rm -rf $dir_name


