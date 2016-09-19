#!/usr/bin/env bash

rm ~/Resources/FunctionalGenomics/data/FANTOM5/processed/CRM/*

perl ~/Resources/FunctionalGenomics/src/FANTOM5dataParser.pl
cd ~/Resources/FunctionalGenomics/data/FANTOM5/processed/CRM/
for i in *.bed
do
	sort-bed $i | uniq > tmp
	mv tmp $i
done
