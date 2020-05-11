#!/bin/bash
i=0
f=0
while [ $f -le 234950 ]
do
	convert -quality 100 test_6/test_$f  test_6/yo_$i.jpg;
	i=$(( i + 1 ))
	f=$(( f + 50 ))
	echo $f
done
# convert -delay 20 -quality 100 yo_*.jpg movie.mpg
