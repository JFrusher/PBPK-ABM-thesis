#!/bin/bash

# Check if a directory was provided, otherwise use the current directory
TARGET_DIR="${1:-.}"

# Check if initial.svg exists to get dimensions
if [ ! -f "$TARGET_DIR/initial.svg" ]; then
    echo "Error: initial.svg not found in $TARGET_DIR. Cannot determine base dimensions."
    exit 1
fi

echo "Processing SVGs in: $TARGET_DIR"

# 1. Get dimensions from initial.svg
# Using 'magick identify' and rounding to even numbers for video/JPEG compatibility
W=$(magick identify -format "%w" "$TARGET_DIR/initial.svg")
H=$(magick identify -format "%h" "$TARGET_DIR/initial.svg")

# Calculate even dimensions (equivalent to your @expr logic)
W1=$(( 2 * (W / 2) ))
H1=$(( 2 * (H / 2) ))
RESIZE_STR="${W1}!x${H1}!"

echo "Target resize dimension: $RESIZE_STR"

# 2. Loop through all s*.svg files in the directory
for svg in "$TARGET_DIR"/s*.svg; do
    # Skip if no files match the pattern
    [ -e "$svg" ] || continue
    
    # Define output filename
    jpg="${svg%.svg}.jpg"
    
    # Only convert if the JPG doesn't already exist (incremental)
    if [ ! -f "$jpg" ]; then
        echo "Converting: $(basename "$svg") -> $(basename "$jpg")"
        magick convert "$svg" -resize "$RESIZE_STR" "$jpg"
    else
        echo "Skipping: $(basename "$jpg") (already exists)"
    fi
done

echo "Done."