#!/bin/bash
i=0
f=0
while [ $f -le 38300 ]
do
	convert -quality 100 test/test_$f  test/yo_$i.jpg;
	i=$(( i + 1 ))
	f=$(( f + 100 ))
	echo $f
done
# convert -delay 20 -quality 100 yo_*.jpg movie.mpg
