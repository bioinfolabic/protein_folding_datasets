#!/bin/sh
# This is a comment!
echo Hello World	# This is a comment, too!

LANG=en_US
min_temp=0.1
temp_step=0.1
max_temp=2.0

pythonFile=tabela_agregacao_v2.py
keyword=bestCoordinates
fileOut=tabela_agregacao

for temp in $(seq $min_temp $temp_step $max_temp);
do
 mypath=Temperatura_$temp
 for item in 1 2 1b 2b;
 do
  python3 $pythonFile $mypath $item\_$keyword $fileOut\_$item
 done
done
