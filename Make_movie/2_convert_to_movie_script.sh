#!/bin/bash
ffmpeg -framerate 100 -i test_6/yo_%00d.jpg -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p test_6.mp4
