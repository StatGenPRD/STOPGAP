#!/usr/bin/env bash

rm ~/Resources/FunctionalGenomics/data/ENCODE/processed/ChIA-PET/*
rm ~/Resources/FunctionalGenomics/data/ENCODE/processed/DHS/*
rm ~/Resources/FunctionalGenomics/data/ENCODE/processed/TFBS/*

perl ~/Resources/FunctionalGenomics/src/ENCODEdataParser.pl

cd ~/Resources/FunctionalGenomics/data/ENCODE/processed/
for i in ChIA-PET/*.bed DHS/*.bed TFBS/*.bed
do
	sort-bed $i | uniq > tmp
	mv tmp $i
done
