#!/bin/bash
ffmpeg -framerate 2 -i test/yo_%00d.jpg -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p test_one.mp4
