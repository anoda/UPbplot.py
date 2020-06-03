#!/bin/sh

# Run all examples

flist=(Chang2006ggg_Fig2.csv Christiansen2009a.csv Christiansen2009b.csv Gonzalez-Leon2011Geosphere_Fig8I.csv Iida2015ia_080122I01.csv Ishihara2015ia_TableS1_68A-64.csv Ito2020TerraNova_FC.csv Iwano2013ia_TableA8.csv Labrado2015cjes.csv Li2016Lithos_Ash1.csv Noda2017bgsj_KT01.csv)

for i in ${flist[@]}; do
	echo $i
	UPbplot.py -i ${i} -c ${i%.csv}.cfg -f > ${i%.csv}.log
done


for j in 1 2; do
	for i in Hoshi2019chemg_Tf${j}*.cfg; do
		echo $i
		UPbplot.py -i Hoshi2019chemg_Tf${j}.csv -c $i -o ${i%.cfg}.pdf -f > ${i%.cfg}.log
	done
done
