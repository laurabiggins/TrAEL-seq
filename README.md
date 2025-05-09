# TrAEL-seq

A collection of scripts for processing TrAEL-seq data

These are used in a nextflow pipeline to process TrAEL-seq data that can be found here: https://github.com/s-andrews/nextflow_pipelines/blob/master/nf_traelseq 

## TrAEL-seq pre-processing - UMI and Poly-T handling (and TrAEL barcodes)

### TrAELseq_preprocessing.py

This is for versions of the TrAEL method that do not incorporate the inline TrAEL barcodes. 

Raw TrAEL-seq FastQ reads are expected to have the following structure:

```
barcode (UMI) (8bp)    //    PolyT     //      Insert
```

**Step 1:**

The script `TrAELseq_preprocessing.py` removes the first 8bp (UMI) of a read and adds the sequence to the end of the readID (separated with a colon, like so:`:UMISEQUENCE`). The quality information is discarded. Empty spaces in the readID are replaced with `_` to preserve the UMI sequence after mapping.

**Step 2:**

After moving the UMI sequences, the script looks for up to 3 T at the start of the sequence, and removes those. Sequences with more than 3 Ts at the 5' end are clipped a maximum of 3 TTT.   
Two output files are produced: T and noT - sequences are assigned depending on the presence of a T at (original) base 9.   
Base 9 is now trimmed whether it is a T or not. We have previously identified that the quality at that position is often low, so that one base is trimmed to aid with mapping in later steps.

This pre-processing script requires Python 3.

In an older version of the script, all FastQ files in a folder were used.
The current version requires filenames to be passed in as arguments to the script. 

```
./TrAELseq_preprocessing.py *fastq.gz
```

###  TrAELseq_preprocessing_UMIplusBarcode.py

This script is used instead of TrAELseq_preprocessing.py for versions of the TrAEL-seq protocol that incorporate inline TrAEL barcodes. 

This is the expected structure of the FastQ files:

```
UMI (8bp)    //   sample-level barcode  4bp   //    PolyT     //      Insert
```

**Step 1:**

This script removes the first 8bp of a read and adds the sequence to the readID. The 4bp sample level barcode is also removed and is written in to the filename. 
The quality information is discarded.

**Step 2:**

After moving the UMI and sample-level barcode sequences, the script looks for up to 3 T at the start of the sequence (original position 13), and removes those. 
Sequences with more than 3 Ts at the 5' end are clipped a maximum of 3 TTT.
Base 13 is trimmed whether it is a T or not. We have previously identified that the quality at that position is often low, so that one base is trimmed to aid with mapping in later steps.

In addition to splitting the sequences by 4bp barcode, there are also T and noT files (as in TrAELseq_preprocessing.py), dependent on the 
presence of a T at position 13.   
This script should produce 20 output fastq files for 1 input fastq file. 
20 = (9 barcodes + 1 unassigned)*2 for T and noT.

Some files may be empty if the barcodes were not all present.
