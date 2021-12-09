# OptimizeGeneNetworks

## Description 

OptimizeGeneNetworks is a script that takes .gmt formatted file and a STRING database file and generates three types of 
files and histograms of network scores per generation. The first of the three files is one that shows summary statistics
of the genetic algorithm used to create highly connected subnetworks of disease associated genes. The genetic algorithm 
also outputs the histograms as a visual summary of each generation. The second file type is a network text file for the top 10 
highest scoring disease associated subnetworks wherein the nodes and edge weights are specified and can be then visualized
in Cytoscape. The third output file is a rewritten .gmt file and the only difference between this one and the original 
is that each gene is given a score.   

## Workflow

### Arguments 

Here is a list of required and optional arguments that can be accessed by 
typing **python "OptimizeGeneNetworks.py" -h**: 

```text
usage: OptimizeGeneNetworks.py [-h] [--gmt GMT] [--sdb SDB] [--ps PS] [--nb NB] [--np NP]

Uses a genetic algorithm to create highly connected disease networks and outputs 
three kinds of files of the results, and displays frequency plots.

optional arguments:
  -h, --help            show this help message and exit
  --gmt GMT, -gmt_formated_file GMT
                        experimental loci file path (default Input.gmt.txt)
  --sdb SDB, -string_database SDB
                        gene interaction file path (default STRING.txt)
  --ps PS, -pop_size PS
                        size of a population of subnetworks
  --nb NB, -n_bins NB   number of bins wanted to create noninformative loci
  --np NP, -n_pops NP   number of populations wanted
```
### Inputs 

Requires one file STRING.txt database, a .gmt format text file, an integer for population size (default 5000), an integer
for the number of bins you want (default 128), and an integer for the number of populations you want (default 1000). 

Input .gmt format file example: 

```text
Fanconi anemia locus 0	Locus for PALB2	NUPR1	CTB-134H23.2	SLC5A11	KIAA0556	CD19	SH2B1	CCDC101	GTF3C1	IL27	ARHGAP17	ERN2	DCTN5	NSMCE1	AQP8	RABEP2	XPO6	ATP2A1	CHP2	BOLA2	KDM8	EIF3C	ATXN2L	LAT	ZKSCAN2	SULT1A1	HS3ST4	EIF3CL	TUFM	NPIPL1	SNX29P2	IL21R	PRKCB	SPNS1	TNRC6A	CACNG3	PLK1	RBBP6	NFATC2IP	APOBR	IL4R	PALB2	SULT1A2	CTD-3203P2.2	GSG1L	SBK1	LCMT1
Fanconi anemia locus 1	Locus for FANCF	CSTF3	FBXO3	SLC17A6	CCDC73	CAPRIN1	RCN1	BDNF	METTL15	CCDC34	EIF3M	LUZP2	BBOX1	CAT	PRRG4	SLC5A12	QSER1	AC103801.2	TCP11L1	SVIP	CD59	NAT10	C11orf91	KCNA4	FIBIN	ARL14EP	ABTB2	LMO2	ELP4	MUC15	DNAJC24	GAS2	LGR4	RP11-17A1.2	MPPED2	KIF18A	FANCF	FSHB	HIPK3	PAX6	DEPDC7	IMMP1L	KIAA1549L	WT1	DCDC5	AC132216.1	ANO5	ANO3	ELF5	EHF	LIN7C
Fanconi anemia locus 2	Locus for RAD51C, BRIP1	CD79B	CACNG1	TANC2	SMG8	RP11-15E18.4	TEX2	YPEL2	RGS9	C17orf72	STRADA	DDX42	TACO1	ICAM2	APOH	PRKCA	FTSJ3	TBX4	DCAF7	GH1	GDPD1	CTD-2510F5.6	METTL2A	MRC2	MAP3K3	PRR11	MED13	C17orf64	TBX2	POLG2	SMURF2	AXIN2	CEP95	INTS2	RAD51C	PPM1E	CA4	CEP112	SMARCD2	C17orf82	USP32	KCNH6	CACNG4	CSH1	RPS6KB1	CTD-2535L24.2	RNFT1	BCAS3	LIMD2	NACA2	RP11-51F16.8	DDX5	APPBP2	SKA2	TRIM37	SCN4A	PTRH2	DHX40	RP11-178C3.1	CCDC47	GNA13	GH2	CSH2	CYB561	HEATR6	VMP1	PSMC5	CSHL1	EFCAB3	TUBD1	CACNG5	BRIP1	PPM1D	AC005544.1	LRRC37A3	CLTC	ERN1	MARCH10	TLK2
Fanconi anemia locus 3	Locus for FANCC	TMOD1	MSANTD3-TMEFF1	TEX10	HIATL2	LPPR1	MRPL50	CDC14B	FOXE1	ZNF510	BAAT	COL15A1	CORO2A	HABP4	INVS	SLC35D2	C9orf156	HSD17B3	ZNF189	GABBR2	C9orf174	FAM22G	ZNF782	ANP32B	XPA	DKFZP434H0512	TSTD2	SEC61B	NR4A3	ANKS6	PTCH1	TMEFF1	MURC	CTSL2	AAED1	STX17	GALNT12	ERP44	TGFBR1	TDRD7	TBC1D2	FANCC	HEMGN	ALG2	MSANTD3	TRIM14	NANS	NCBP1	KRT8P11	ERCC6L2	ZNF367
Fanconi anemia locus 4	Locus for FANCA	DBNDD1	AC133919.6	RP11-356C4.2	MC1R	SPIRE2	C16orf3	CENPBD1	FANCA	DEF8	GAS8	PRDM7	TCF25
...
```
Input STRING.txt database example: 
```
ARF5	DVL2	0.166000
ARF5	DYRK4	0.166000
ARF5	PPP5C	0.254968
ARF5	MAP4K5	0.157276
ARF5	RALBP1	0.156000
ARF5	PKP2	0.160210
ARF5	ACAP1	0.328000
ARF5	MAP2K5	0.242000
ARF5	MYO15A	0.272395
ARF5	MAPK13	0.190000
ARF5	STX1B	0.263160
ARF5	MAPK12	0.190000
ARF5	MAPK1	0.190000
ARF5	MYH9	0.252822
ARF5	SOS2	0.199000
ARF5	EIF5	0.214358
ARF5	PABPC1L	0.196000
ARF5	CORO1A	0.163000
ARF5	PARD6A	0.194000
ARF5	PDIA2	0.157397
...
```
### Command

To run this program in the command line interface type: 
```text
python "OptimizeGeneNetworks.py" --gmt "Input.gmt.txt" --sdb "STRING.txt" --ps 5000 --nb 128 --np 1000
```
You can replace --gmt and --sdb with your own desired input files as well as change 
the population size, the number of bins and the number of populations. The different output files will automatically 
appear in the current working directory. 

To run this program interactively type: 

```python
OptimizeGeneNetworks("Input.gmt.txt", "STRING.txt", 5000, 128, 1000)
```
If the input files are located in a different directory then you can put their respective file paths. 

### Output 

Example of Genetic Algorithm summary file: 

```text
Generation: 0
Selection score's mean: 0.225712
Selection score's standard deviation: 0.49267210009867934
Selection score's variance: 0.24272579821564313

Generation: 1
Selection score's mean: 0.520972
Selection score's standard deviation: 0.793680969352436
Selection score's variance: 0.6299294811122225

Generation: 2
Selection score's mean: 0.767872
Selection score's standard deviation: 1.0413864646564954
Selection score's variance: 1.084485768769754
...
```

Example of network text file formatting:

```text
PLK1	RP11-298P3.4	0.807000
PLK1	ERCC4	0.434883
PLK1	BLM	0.297000
FANCF	FANCC	0.984000
FANCF	FANCA	0.999000
FANCF	FANCE	0.982000
FANCC	FANCF	0.984000
FANCC	FANCA	0.999000
FANCC	FANCD2	0.807000
FANCC	FANCE	0.999000
FANCA	FANCF	0.999000
FANCA	FANCC	0.999000
FANCA	FANCE	0.999000
FANCA	ERCC4	0.566000
FANCA	BLM	0.906000
PIK3C2B	ERCC4	0.248000
PIK3C2B	BLM	0.304000
FANCD2	FANCC	0.807000
...
```

Example of output .gmt file:

```text
Fanconi anemia locus 0	Locus for PALB2	NUPR1 124.74	CTB-134H23.2 NA	SLC5A11 3.68	KIAA0556 0	CD19 19.62	SH2B1 0	CCDC101 260.83	GTF3C1 0	IL27 0	ARHGAP17 0	ERN2 341.64	DCTN5 0	NSMCE1 588.62	AQP8 24.08	RABEP2 6.22	XPO6 178.03	ATP2A1 430.22	CHP2 491.35	BOLA2 271.7	KDM8 29.71	EIF3C 929.46	ATXN2L 159.81	LAT 0	ZKSCAN2 93.56	SULT1A1 24.47	HS3ST4 0	EIF3CL 563.54	TUFM 1421.87	NPIPL1 NA	SNX29P2 NA	IL21R 0	PRKCB 668.16	SPNS1 57.48	TNRC6A 0	CACNG3 0	PLK1 5533.09	RBBP6 73.65	NFATC2IP 2872.36	APOBR NA	IL4R 29.38	PALB2 711.64	SULT1A2 0	CTD-3203P2.2 NA	GSG1L 0	SBK1 92.38	LCMT1 58.17
Fanconi anemia locus 1	Locus for FANCF	CSTF3 974.77	FBXO3 0	SLC17A6 19.12	CCDC73 NA	CAPRIN1 18.33	RCN1 91.66	BDNF 0	METTL15 163.89	CCDC34 0	EIF3M 714.2	LUZP2 0	BBOX1 48.0	CAT 196.3	PRRG4 0	SLC5A12 15.64	QSER1 0	AC103801.2 NA	TCP11L1 30.08	SVIP 0	CD59 0	NAT10 656.18	C11orf91 NA	KCNA4 28.72	FIBIN NA	ARL14EP 65.38	ABTB2 0	LMO2 138.47	ELP4 39.03	MUC15 NA	DNAJC24 201.47	GAS2 0	LGR4 452.45	RP11-17A1.2 NA	MPPED2 0	KIF18A 399.5	FANCF 18129.54	FSHB 0	HIPK3 712.36	PAX6 45.25	DEPDC7 100.67	IMMP1L 20.21	KIAA1549L 0	WT1 261.77	DCDC5 NA	AC132216.1 NA	ANO5 0	ANO3 NA	ELF5 192.85	EHF 40.2	LIN7C 8.85
Fanconi anemia locus 2	Locus for RAD51C, BRIP1	CD79B 37.41	CACNG1 0	TANC2 88.4	SMG8 174.94	RP11-15E18.4 NA	TEX2 5.78	YPEL2 18.78	RGS9 0	C17orf72 NA	STRADA 0	DDX42 44.04	TACO1 58.13	ICAM2 0	APOH 10.05	PRKCA 688.73	FTSJ3 396.25	TBX4 26.9	DCAF7 47.6	GH1 24.34	GDPD1 60.22	CTD-2510F5.6 NA	METTL2A 36.67	MRC2 0	MAP3K3 825.56	PRR11 203.2	MED13 15.47	C17orf64 1841.01	TBX2 9.39	POLG2 870.17	SMURF2 193.49	AXIN2 71.97	CEP95 23.45	INTS2 0	RAD51C 3980.63	PPM1E 445.11	CA4 0	CEP112 0	SMARCD2 204.61	C17orf82 0	USP32 55.24	KCNH6 13.28	CACNG4 NA	CSH1 19.58	RPS6KB1 1249.4	CTD-2535L24.2 NA	RNFT1 46.57	BCAS3 0	LIMD2 0	NACA2 1044.59	RP11-51F16.8 0	DDX5 665.48	APPBP2 228.04	SKA2 14.54	TRIM37 47.09	SCN4A 0	PTRH2 71.98	DHX40 75.75	RP11-178C3.1 696.21	CCDC47 98.55	GNA13 411.21	GH2 0	CSH2 19.58	CYB561 0	HEATR6 0	VMP1 0	PSMC5 529.32	CSHL1 0	EFCAB3 120.36	TUBD1 406.83	CACNG5 0	BRIP1 1471.02	PPM1D 448.81	AC005544.1 NA	LRRC37A3 0	CLTC 291.09	ERN1 212.17	MARCH10 126.79	TLK2 9.01
Fanconi anemia locus 3	Locus for FANCC	TMOD1 22.72	MSANTD3-TMEFF1 NA	TEX10 199.87	HIATL2 18.56	LPPR1 0	MRPL50 0	CDC14B 1385.5	FOXE1 231.68	ZNF510 114.11	BAAT 0	COL15A1 0	CORO2A 250.56	HABP4 240.29	INVS 88.4	SLC35D2 42.93	C9orf156 0	HSD17B3 0	ZNF189 114.25	GABBR2 0	C9orf174 0	FAM22G NA	ZNF782 114.11	ANP32B 230.84	XPA 2380.11	DKFZP434H0512 NA	TSTD2 14.02	SEC61B 89.22	NR4A3 109.93	ANKS6 288.76	PTCH1 26.0	TMEFF1 NA	MURC 0	CTSL2 20.48	AAED1 0	STX17 268.53	GALNT12 0	ERP44 254.49	TGFBR1 380.1	TDRD7 0	TBC1D2 432.37	FANCC 20379.6	HEMGN 0	ALG2 212.47	MSANTD3 0	TRIM14 0	NANS 14.93	NCBP1 1179.33	KRT8P11 4.52	ERCC6L2 122.05	ZNF367 114.11
Fanconi anemia locus 4	Locus for FANCA	DBNDD1 NA	AC133919.6 0	RP11-356C4.2 NA	MC1R 124.43	SPIRE2 0	C16orf3 NA	CENPBD1 0	FANCA 27121.77	DEF8 25.62	GAS8 0	PRDM7 0	TCF25 14.31
...
```

Example of genetic algorithm histogram: 

![Image of histogram example](https://github.com/acolorado1/Module4Day3/blob/1e1a56674ce582d4363ac421c76a7e70b3b92dcf/Module4Day3_hist_example.png)

Note: If script is run through the command line these histograms will not appear. 

### Visualization 
Visualization was done in Cytoscape. Note, when importing the network into Cytoscape as it is a text file one will have
to click Selct Network, then click column 1 and select source node, then click column 2 and select target node, and 
click column 3 and select edge attribute. 

The following is the visualization of the example network text file: 

![Image of subnetwork example](https://github.com/acolorado1/Module4Day3/blob/545e3f6db9c289c043c2c1a3f4e3134df4aaa064/Day3_Output_Network1_pvalue0.0.PNG)

## Installation and Dependencies
You must have Python 3 installed. Any Python 3 version should work but it was written in Python 3.9 using a Windows-based 
operating system. Packages random, argparse 1.4.0, the functions mean, stdev, and variance from the statistics package, 
and matplotlib 3.5.0 will need to be installed. 

## Contact 
Angela Sofia Burkhart Colorado - angelasofia.burkhartcolorado@cuanschutz.edu

Project Link: https://github.com/acolorado1/Module4Day3.git