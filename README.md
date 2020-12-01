# Incubation Predictor and Mutation Rates Calculator
The code provided in this repository allows the user to calculate mutation rates and predict ssRNA viruses incubation times, as described in [1]

### Citation
[1] Ayal B. Gussow*, Noam Auslander*#, Yuri I. Wolf, Eugene V. Koonin# [Prediction of the incubation period for COVID-19 and future virus disease outbreaks](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-020-00919-9) BMC Biology, 2020
(*) These authors contributed equally, (#) Corresponding authors


### Mutation rate calculation
calc_mutation_rate.py receives two file paths as input: 

1. Path to Fasta file of alignment.
2. Path to tree file in Newick format. The tree must be based on the alignment file, thus the organism names should 
match between the two files. 

### Incubation Prediction
calc_incubation.py receives four input variables:
 
1) Path to Fasta of aligned sequence. 
2) Path to GenBank file.
4) The reference genome id. 
5) An optional argument 'training_families', describing which model to use based on the families used in training. 
If not provided, default is Coronaviridae. Other options are Coronaviridae_Paramyxoviridae and Coronaviridae_Pneumoviridae 
(training on <em>Coronaviridae</em> plus either <em>Paramyxoviridae</em> or <em>Pneumoviridae</em>, respectively). 

The reference genome must appear in the Fasta file and GenBank file.

The scripts can be run from the command line or loaded into a Python script.

The scripts were tested on Python 3.7.1 and 3.7.4, with the package versions in requirements.txt.

## Importing Into a Python Script

Below are examples of how to load the scripts as a Python module (needs to be run from within the same directory as the incubation Python scripts).

### Mutation rate calculation
To incorporate into a Python script, use the calc_mutation_rate function.
For example, using the sample data in this repository:

```python
import calc_mutation_rate

# setup the data and calculate mutation rate
mutation_rate, mutation_counts, merged_array = calc_mutation_rate.calc_mutation_rate(
    "sample_data/229e.afa",
    "sample_data/229e.tre"
)

# print the mutation rate
print("The mutation rate is: {}".format(mutation_rate))
```

### Incubation Prediction
To incorporate into python script, use the calc_incubation function. 
For example, using the sample data in this repository:

```python
import calc_incubation

# setup the model data and run the prediction
pred_incubation = calc_incubation.calc_incubation(
    "sample_data/229e.afa",
    "sample_data/229e.gb",
    "NC_002645.1"
)

# print the prediction
print("The predicted incubation is: {}".format(round(pred_incubation, 3)))
```

## Command Line Usage
Below see sample commands to run the scripts using the data files in the sample\_data directory.

### Mutation Rate Calculation
To calculate mutation rate:

```python calc_mutation_rate.py --alignment sample_data/229e.afa --tree sample_data/229e.tre```

### Incubation prediction
To calculate incubation rate with the default model trained on <em>Coronaviridae</em>:

```python calc_incubation.py --alignment sample_data/229e.afa --genbank sample_data/229e.gb  --genome-id NC_002645.1```

Alternatively, to calculate incubation rate with the model trained on <em>Coronaviridae</em> and <em>Pneumoviridae</em>:

```python calc_incubation.py --alignment sample_data/229e.afa --genbank sample_data/229e.gb --genome-id NC_002645.1 --training_families Coronaviridae_Pneumoviridae```


