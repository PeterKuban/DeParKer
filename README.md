# Distributed Parallel K-mers Counter: DeParKer

Copyright 2016 Peter Kubán, Mária Lucká and Tomaáš Farkaš <br />
peter_kuban (at) stuba.sk

## DeParKer
DeParKer is a software for counting _k_-mers usable on parallel distributed systems (using OpenMPI and OpenMP).

## Compilation
To compile the source code, type

```
make
``` 

It will create two executable files in bin directory:

* bin\deparker generates binary files containing _k_-mers with their counts, 
* bin\dump_counted produces human readable (tab delimited) file having _k_-mers with their counts.

## Usage
To go to bin directory:

	cd bin

Run deparker at first, and then run dump_counted. dump_counted uses deparker generated files as input, and outputs list of human readable _k_-mers with their counts.

To run deparker, use

```
Usage: mpirun -n <num_of_processors> deparker [options] <input_file_path>
-	input file can be in FASTA format or as raw strings (reads)
Options (default value in (), * - required):
	-k=uint32	*Length of k-mer (31)
	-t=uint32	Number of threads (2)
	-m=uint32	Available memory for core/processor in MB (2048)
	-o=string	*Output file name (full path)
	-l=uint32	Don't output k-mers with count lower than this value (2)
```
Example: ```mpirun -n 4 deparker -t 2 -k 31 -l 5 -o /path/to/output/ input_file.fa```

To run dump_counted, use
```
Usage: dump_kmers [options] <input_file_path>
Options (* - required):
	-o=string	*Output file name (full path)

```
Example: ```dump_counted -o /path/to/output/output.txt kmers_counted_0_1234```


## Contact
For questions, suggestions, bugs, or other issues, please contact:

```
Peter Kubán
peter_kuban (at) stuba.sk
```
