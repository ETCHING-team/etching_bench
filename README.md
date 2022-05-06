# A benchmark tool for ETCHING

### Version 0.2.12

This program is a benchmarking tool implementing SURVIVOR for somatic structural variation (SV) callers. 

##### Supporting SV callers

* ETCHING, DELLY, LUMPY, Manta, SvABA, novoBreak, and GRIDSS

**Note**: This is not a general-purpose tool. We did not check other callers.



## Installation

##### Requirement

* cmake >=3.11 (We did not check <3.11)

##### installation

```bash
mkdir build
cd build
cmake ../
make
```

Then, you can find  ```etching_bench```  in the path ```build```.

##### Test

```bash
cd ../test/
../build/etching_bench -c test.conf
```



## Usage

```bash
Usage: etching_bench [options]

Required:
	-c FILE    Config file (required)

Options:
	-o STR     Outfile prefix [etching_bench]
	-s INT     Consensus cutoff for silver standard [3]
	           TRUE if detected by >=3 (in default) callers.
	-t FILE    Truth set
	-w INT     Merge window size [10]
	-m INT     Minimun SV size [100]
	-M INT     Maximum SV size
	-I         Remove IMPRECISE SVs for all callers
	-X         Exclude own prediction in a silver standard set [default]
	-K         Keep own prediction in a silver standard set
	-F         Calculate performance at F1-score-maximizing cutoff
	-R FILE    Re-use a given.benchmark.annotated.vcf (skip to calculation)
	--intra    Only intra-chromosomal SVs (DEL, DUP, or INV)
	--inter    Only inter-chromosomal SVs (TRA)
	--version  Print version
	--version  Print version
	-h         Print this message

Note: Not for general use yet.
      We guarantee only for the programs listed below:
      ETCHING, DELLY, LUMPY, Manta, SvABA, novoBreak, and GRIDSS

```



#### 1. Edit a ```conf``` file. 

```bash
vi input.conf
```

The format of the ```conf``` file

```bash
TOOL_NAME	VCF_FILE_NAME	[cutoff]	["noimprecise"]
```

The third and fourth columns are optional, and their positions are interchangeable. 

An example of the ```conf``` file:

```
etching	ETCHING.vcf	0.4
delly	DELLY.vcf	noimprecise	18
Lumpy	LUMPY.vcf	12	noimprecise
Manta	manta.vcf	noimprecise
SvABA	svaba.vcf	0
novobreak	novoBreak.vcf	40
GRIDSS	Gridss.vcf	400
```

**Note**: This program converts all letters in TOOL_NAME to lower case. Thus, ETCHING, Etching, and even eTChiNg are acceptable.



##### Supporting tools

| Caller    | Default cut-off | Score to be used                        |
| --------- | --------------- | --------------------------------------- |
| ETCHING   | 0.4             | QUAL field value                        |
| DELLY     | 18              | DV + RV in tumor sample of FORMAT field |
| LUMPY     | 12              | SU in tumor sample of FORMAT field      |
| Manta     | 40              | SOMATICSCORE in INFO field              |
| SvABA     | 0               | QUAL field value                        |
| novoBreak | 40              | QUAL field value                        |
| GRIDSS    | 400             | QUAL field value                        |



##### About the option, "noimprecise" for each tool

This option excludes ```IMPRECISE``` SVs for the tool in benchmarking. If you want to exclude all ```IMPRECISE``` SVs for all callers, use the ```-I``` option instead in running ```etching_bench```. 

**Note**: Removing ```IMPRECISE``` will increase precision, while it may drop recall (or sensitivity).



#### 2. Run etching_bench 

```bash
etching_bench -c input.conf [-o output_prefix] [options]
```

##### Available options

	-s INT     Consensus cutoff for silver standard [3]
	           TRUE if detected by >=3 (in default) callers.
	-w INT     Merge window size [10]
	-m INT     Minimun SV size [100]
	-M INT     Maximum SV size
	-I         Remove IMPRECISE SVs for all callers
	-X         Exclude own prediction in a silver standard set [default]
	-K         Keep own prediction in a silver standard set
	-F         Calculate performance at F1-score-maximizing cutoff
	-R FILE    Re-use a given.benchmark.annotated.vcf (skip to calculation)
	--intra    Only intra-chromosomal SVs (DEL, DUP, or INV)
	--inter    Only inter-chromosomal SVs (TRA)



##### About silver-standard

Our silver-standard calls the SVs predicted by N (3 in default) *other callers* as TRUE. It means ETCHING's silver-standard has no ETCHING, and there is no DELLY in DELLY's silver-standard, and so on. If you include the target tool (using ```-K```), it may cause a bias toward a tool of very high sensitivity with very low precision. Thus, we do not recommend it but keep the default (```-X``` to specify it).



#### 3. Output files

| Suffix of output file         | Description                                            |
| ----------------------------- | ------------------------------------------------------ |
| benchmark.txt                 | Final table of performances relevant to their cut-offs |
| TOOL_NAME.performance.ALL.txt | Performances of all SVs for different cut-offs         |
| TOOL_NAME.performance.DEL.txt | Performances of DELs for different cut-offs            |
| TOOL_NAME.performance.DUP.txt | Performances of DUPs for different cut-offs            |
| TOOL_NAME.performance.INV.txt | Performances of INVs for different cut-offs            |
| TOOL_NAME.performance.TRA.txt | Performances of TRAs for different cut-offs            |

You can draw precision-recall (PR) curves using the ```performance.XXX.txt``` and find the area under the PR curve (auPR) at the tail of each ```performance.XXX.txt```.



## Contact

Jang-il Sohn (sohnjangil@gmail.com)
