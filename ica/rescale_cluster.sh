#!/bin/sh


cat cluster_list.txt|while read l
do
clu_num=`echo $l|awk '{print $2}'`
freq=`echo $l|awk '{print $1}'`
echo $clu_num $freq
../rescale.py mask.fits res_${freq}_${clu_num}.fits model_${freq}.fits cluster_${freq}.fits 
done
