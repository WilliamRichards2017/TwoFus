# TwoFus
**Author** - Will Richards

Structural variant detection plugin for  kmer based denovo variant detection software [rufus](https://github.com/jandrewrfarrell/RUFUS). 
Used to call translocations, inversions, mobile elements, and multi-megabase insertions


#  Run with rufus (recommended)

 This plugin uses intermediate files generated from rufus as input, and writes out called variants to the final vcf file produced by rufus.  Rufalu is installed as an external dependency inside rufus, and the best way to use this tool to call alus is to run rufus directly, by follow the install and running instructions [here](https://github.com/jandrewrfarrell/RUFUS)
# Running as standalone app

This plugin is primarily ment to be run inside of rufus, although it can be run as a standalone app if a user has the appropriate files.

### Dependencies
[Cmake](https://cmake.org/download/) - c++ build management

### Install
```
cd bin
cmake ../
make
```

### Run 

To run TwoFus, you must run the aluDetect script with the following positional paramers:
```
./tools/runTwofus $1={{path-to-bam-file}} $2={{path-to-mobile-element-fasta-file}} $3={{path-to-reference-file}} $4={{path-to-jhash-file}} $5+={{paths-to-parent-bam-files}}
```

