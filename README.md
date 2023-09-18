# AlphaFold unmasked <br><sup>to integrate experiments and predictions</sup>

![H1111](fig/header.png)

## Installing AF_unmasked

The installation and setup procedure is the same as for the regular version of AlphaFold (non-docker version). We recommend Anaconda and mamba along with pip3 to manage the necessary software packages:

1. [Install Anaconda/Miniconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

2. Set up conda environment, install dependencies:

```bash
conda create -n AF_unmasked -c conda-forge python=3.9 mamba
conda activate AF_unmasked

# clone this repository
git clone https://github.com/clami66/AF_unmasked.git
cd AF_unmasked/

# install requirements
mamba env update --file environment.yaml
python -m pip install -r requirements.txt
```

3. [optional] Download and set up the AF parameters and sequence databases. We recommend downloading the reduced set of databases since evolutionary inputs are not as important when a good template is provided. If the full databases are needed, run the following by omitting `reduced_dbs`:

```bash
cd scripts
chmod +x download_all_data.sh

./download_all_data.sh ../AF_data/ reduced_dbs
```

If you have databases and parameters from a precedent AlphaFold installation, it is not necessary to repeat this step, just make sure that the paths inside `databases.flag` point to the right directories.

4. [optional] Install [lDDT_align](https://github.com/clami66/lDDT_align) if you want to perform superposition-free structural alignments

## Preparing multimeric templates

This version of AlphaFold comes with a python script `prepare_templates.py` to set up multimeric templates before a run.

**Quick start**

If you have a `.fasta` file containing multiple sequences:

```
>H1137,subunit1|
MTEPPAPTAPLNKPKTPPYKLAGLILGLVGVLVLALTWMQFRGQFEDKVQLTVLSGRAG...
>H1137,subunit2|
MSIKGTLFKLGIFSLVLLTFTALIFVVFGQIRFNRTTEYSAIFKNVSGLRDGQFVRAAG...
>H1137,subunit3|
MRTLQGSDRFRKGLMGVIVVALIIGVGSTLTSVPMLFAVPTYYGQFADTGGLNIGDKVR...
...
```

And a `.pdb`/`.cif` file containing as many template chains as there are target chains in the fasta file, you can run e.g.:

```
python prepare_templates.py --target examples/H1137/H1137.fasta \
    --template examples/H1137/H1137.pdb \
    --output_dir AF_models/ \
    --align
```

When using a `.fasta` target, the outputs will be saved in a subfolder inside `output_dir` with the same name as the fasta file. In this case, the outputs will be stored inside `AF_models/H1137` because the fasta filename is `H1137.fasta`.

**Chain mapping flags**

The previous example assumes that the first chain in the `.fasta` file maps to the first chain in the `.pdb` template, and so on. If not, it is necessary to specify the mapping. For example, if the first sequence in the `.fasta` file maps to the `B` chain in the template, and the second sequence maps to the `C` chain in the template, you can run:

```
python prepare_templates.py --target examples/H1142/H1142.fasta \
    --template examples/H1142/H1142.pdb \
    --output_dir AF_models/ \
    --align \
    --target_chains A B \
    --template_chains B C
```

The `--target_chains`/`--template_chains` mapping flags are also necessary when the template contains more/fewer chains than there are sequences in the input `.fasta` file.

**Input structures instead of sequences**

It is possible to prepare templates starting from `.pdb`/`.cif` files instead of `.fasta` sequences. This is useful, for example, when the user wants to start from monomeric predictions of each target chain and align them against a multimer template structure. For example, if two unbound structures for chains `A` and `B` are available in PDB format:

```
python prepare_templates.py --target examples/H1142/casp15_predictions/unbound_chain_A.pdb \
    examples/H1142/casp15_predictions/unbound_chain_B.pdb \
    --template examples/H1142/H1142.pdb \
    --output_dir AF_models/H1142 \
    --align
```

NB: in this case, the user needs to manually add to the output directory the name of the `.fasta` file that will be used in the AlphaFold run (`H1142.fasta` -> `--output_dir AF_models/H1142`).

The target chains can also come from the same PDB file, in that case it might be necessary to provide the chain mapping flags:

```
python prepare_templates.py --target examples/H1142/casp15_predictions/unbound_chains.pdb \
    --template examples/H1142/H1142.pdb \
    --output_dir AF_models/H1142 \
    --align \
    --target_chains A B \
    --template_chains B C
```

The target/template files are in `.pdb` format by default, but mmCIF is also supported. The `--mmcif_target`/`--mmcif_template` flags are expected in that case.

**Assembling unbound structures onto template**

When a template for an interaction is available that is a remote homolog of the target interaction, it might be useful to superimpose unbound monomers onto the template to create a coarse interaction model. Then, the coarse model made from putting together the unbound monomers will be itself used as template. This can be done with the `--superimpose` flag:

```
python prepare_templates.py --target examples/H1142/casp15_predictions/unbound_chain_A.pdb \
    examples/H1142/casp15_predictions/unbound_chain_B.pdb \
    --template examples/H1142/H1142.pdb \
    --output_dir AF_models/H1142 \
    --align \
    --superimpose
```

**Adding further templates**
AlphaFold takes up to four structural templates as input. Once the first template has been generated with `prepare_templates.py`, three more can be added with the `--append` flag. These can be the same template as the first, repeated three more times, or different templates:

```
# prepare the first template
python prepare_templates.py --target examples/H1137/H1137.fasta --template examples/H1137/H1137.pdb --output_dir AF_models/ --align
# running the same command three more times to fill the four template slots while using different templates
python prepare_templates.py --target examples/H1137/H1137.fasta --template examples/H1137/H1137.pdb --output_dir AF_models/ --align --append
python prepare_templates.py --target examples/H1137/H1137.fasta --template examples/H1137/H1137.pdb --output_dir AF_models/ --align --append
python prepare_templates.py --target examples/H1137/H1137.fasta --template examples/H1137/H1137.pdb --output_dir AF_models/ --align --append
```

## Outputs

Having run `prepare_templates.py` four times (one per template), the output directory `AF_models/` will look as follows:

```bash
AF_models/H1137/
├── H1137.fasta  # target fasta file
├── H1137.pdb    # template pdb file
├── msas
│   ├── A
│   │   └── pdb_hits.sto
│   ├── B
│   │   └── pdb_hits.sto
...
└── template_data
    ├── mmcif_files
    │   ├── 0000.cif
    │   ├── 0001.cif
    │   ├── 0002.cif
    │   └── 0003.cif
    ├── pdb_seqres.txt
    └── templates.flag
```

The `template_data/` subfolder mimics AlphaFold's database of PDB structures, where only four PDB structures are included (one per template). Whatever the name of the template, the `.cif` PDB files are renumbered from 0000 to 0003. 

If the `--align` option has been used, then the `msas/` folder will contain an alignment file for each chain matched in the template (`pdb_hits.sto`). These should replace the ones generated by a regular AlphaFold run. If the `--align` option has not been used, then AlphaFold will search for templates in the `template_data/pdb_seqres.txt` sequence database and generate the `pdb_hits.sto` files on its own.

The `template_data/templates.flag` file is a flagfile that should be passed to AlphaFold when it's time to perform a prediction. It cointains the information AlphaFold needs to find all the necessary template/alignment information as it's been generate by `prepare_templates.py`.

## Running AlphaFold

Once templates have been prepared, invoke AlphaFold with the generated flagfile (inside the `template_data` folder) along with the standard flagfile (`databases.flag` in this repository).

Use the `--cross_chain_templates` or `--cross_chain_templates_only` flags if you want to use both intra- and inter-chain constraints from the template, or inter-chain constraints alone:

```
python run_alphafold.py --fasta_paths examples/H1137/H1137.fasta \
    --flagfile ./databases.flag \
    --flagfile examples/H1137/template_data/templates.flag \
    --output_dir AF_models \
    --cross_chain_templates \
    --dropout \
    --model_preset='multimer_v2'
```

## Predicting homomers

whenever running with homomers, or multimers containing multiple copies of any given chain, make sure to add the `--separate_homomer_msas` flag, in order to force AlphaFold to read the correct `pdb_hits.sto` template alignment:

```
python run_alphafold.py --fasta_paths examples/H1137/H1137.fasta \
    --flagfile ./databases.flag \
    --flagfile examples/H1137/template_data/templates.flag \
    --output_dir AF_models \
    --cross_chain_templates \
    --dropout \
    --model_preset='multimer_v2' \
    --separate_homomer_msas
```

## Clipping the MSAs to speedup computation

Use the `[uniprot,mgnify,uniref,bfd]_max_hits` flags to limit the number of sequences to include from each alignment file. For example, if we only want to use 200 sequences from mgnify and uniref, while only keeping a single sequence from other alignments:

```
python run_alphafold.py --fasta_paths examples/H1137/H1137.fasta \
    --flagfile ./databases.flag \
    --flagfile examples/H1137/template_data/templates.flag \
    --output_dir AF_models \
    --cross_chain_templates \
    --dropout \
    --model_preset='multimer_v2' \
    --separate_homomer_msas \
    --uniprot_max_hits 1 \
    --mgnify_max_hits 200 \
    --uniref_max_hits 200 \
    --bfd_max_hits 1
```

Limiting the number of alignments from MSAs forces AlphaFold to rely more on the templates, while speeding up computation. Maximum speedup is achieved by making sure that the maximum total number of sequences in the final alignment is no more than 512.

## References

If you use AF_unmasked you can reference: [DOI]

As well as the original [AlphaFold](https://doi.org/10.1038/s41586-021-03819-2) and [AlphaFold-Multimer](https://www.biorxiv.org/content/10.1101/2021.10.04.463034v1) papers.
