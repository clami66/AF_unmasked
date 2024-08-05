import os
import re
import sys
import pickle
import tempfile
import traceback
import subprocess
from pathlib import Path
from string import ascii_uppercase, ascii_lowercase
from argparse import ArgumentParser
import numpy as np
from Bio import Align, SeqIO
from Bio.Data import IUPACData
from Bio.SeqUtils import seq1, seq3
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.kdtrees import KDTree
from Bio.PDB import MMCIFParser, PDBParser, Superimposer

ascii_upperlower = ascii_uppercase + ascii_lowercase


def parse_args():
    parser = ArgumentParser(
        description="Align and format templates of interest to your target protein structures. The outputs can be used to run AlphaFold on the selected templates"
    )
    parser.add_argument(
        "--target",
        metavar="targets",
        type=str,
        nargs="+",
        help="path to target(s) files. Multiple target files can be passed, e.g. one per chain in the template. Make sure you specify which chains should be aligned with --target_chains.",
        required=True,
    )
    parser.add_argument(
        "--template",
        metavar="template",
        type=str,
        help="path to template file. Only one template file should be passed",
        required=True,
    )
    parser.add_argument(
        "--out_dir", "--output_dir",
        type=str,
        help="path to output folder, which will include the template's structure, stockholm alignments and flagfile for usage on AlphaFold",
        default="./AF_models",
    )
    parser.add_argument(
        "--align",
        action="store_true",
        default=True,
        help="whether the template alignment file (pdb_hits.sto) should be generated",
    )
    parser.add_argument(
        "--noalign",
        action="store_false",
        dest="align",
        help="turn alignments off, let AlphaFold generate the pdb_hits.sto files",
    )
    parser.add_argument(
        "--align_tool",
        type=str,
        help="alignment tool",
        choices=["blast", "lddt_align", "tmalign"],
    )
    parser.add_argument(
        "--revision_date",
        type=str,
        help="revision date to set on output mmcif",
        default="2100-01-01",
    )
    parser.add_argument(
        "--mmcif_target",
        default=False,
        action="store_true",
        help="The target(s) are in mmCIF format",
    )
    parser.add_argument(
        "--mmcif_template",
        default=False,
        action="store_true",
        help="The template is in mmCIF format",
    )
    parser.add_argument(
        "--append",
        default=False,
        action="store_true",
        help="Add template alignment to existing 'pdb_hits.sto' alignments",
    )
    parser.add_argument(
        "--superimpose",
        default=False,
        action="store_true",
        help="Superimpose target chains to template and use the superposition as template structure instead",
    )
    parser.add_argument(
        "--inpaint_clashes",
        default=True,
        action="store_true",
        help="If clashing residues between chains should be deleted from the template, this will let AF inpaint them during the modelling stage",
    )
    parser.add_argument(
        "--noinpaint_clashes",
        dest="inpaint_clashes",
        action="store_false",
        help="Turn automatic inpainting of clashes off",
    )
    parser.add_argument(
        "--target_chains",
        metavar="target_chain",
        type=str,
        nargs="+",
        help="pdb chains for target",
    )
    parser.add_argument(
        "--template_chains",
        metavar="target_chain",
        type=str,
        nargs="+",
        help="pdb chains for template",
    )

    return parser.parse_args()


def is_fasta(path):
    records = list(SeqIO.parse(path, "fasta"))
    return True if records else False


def load_PDB(path, n_model=0, is_mmcif=False):

    if not is_mmcif:
        pdb_parser = PDBParser(QUIET=True)
    else:
        pdb_parser = MMCIFParser(QUIET=True)

    try:
        structure = pdb_parser.get_structure("-", path)
        model = structure[n_model]
    except Exception as e:
        print("ERROR: is the file in the correct format? (.pdb, .cif)")
        if not is_mmcif:
            print("       (use -mmcif_model or -mmcif_native with mmCIF inputs)")
        print(traceback.format_exc())
        sys.exit(1)
    return model


def remove_extra_chains(model, chains_to_keep):
    chains = [chain.id for chain in model.get_chains()]

    chains_to_remove = set(chains).difference(set(chains_to_keep))
    for chain in chains_to_remove:
        model.detach_child(chain)


def remove_hetatms(model):
    chains = [chain.id for chain in model.get_chains()]
    residues_to_delete = []

    for chain in chains:
        residues = model[chain].get_residues()

        for res in residues:
            if res.id[0] != " ":
                residues_to_delete.append(res.get_full_id())
    for _, _, chain, res in residues_to_delete:
        model[chain].detach_child(res)


def detect_and_remove_clashes(model, clash_threshold=3.5):

    # Extract CA atoms from the structure
    ca_atoms = []
    ca_data = []
    for chain in model:
        for residue in chain:
            if residue.has_id('CA'):
                coords = residue['CA'].get_coord()
                ca_atoms.append(coords)
                ca_data.append((chain, int(residue.id[1])))

    cb_tree = KDTree(np.array(ca_atoms, dtype="d"))
    # Check for clashes and delete clashing residues
    to_delete = set()
    for i, coord in enumerate(ca_atoms):
        clashing_indices = cb_tree.search(np.array(coord, dtype="d"), clash_threshold)
        for j in clashing_indices:
            if ca_data[i][1] != ca_data[j.index][1]:
                to_delete.add(ca_data[i])
                to_delete.add(ca_data[j.index])

    n_deleted = 0
    for chain in model:
        for i in to_delete:
            if i[0] == chain:
                try:
                    chain.detach_child((' ', i[1], ' '))
                    n_deleted += 1
                except KeyError:
                    pass
    return n_deleted, model

def get_fastaseq(model, chain):
    return "".join(seq1(aa.get_resname()) for aa in model[chain].get_residues())


def write_seqres(path, sequences, chains, seq_id="0000", append=False):
    """
    Format:
    >0000_A mol:protein length:216
    QSALTQPASVSGSPGQSITISCTGTSSDVGGY ...
    >0000_B ...
    ...
    """
    with open(path, mode="a" if append else "w") as out:
        for sequence, chain in zip(sequences, chains):
            out.write(f">{seq_id}_{chain} mol:protein length:{len(sequence)}\n")
            out.write(f"{sequence}\n")
    return


def format_alignment_stockholm(alignments, hit_id="", hit_chain=""):
    ref_align, query_align = alignments
    len_ref_align = len([aa for aa in ref_align if aa != "-"])
    # GS record, header
    formatted_alignment = f"#=GS {hit_id}_{hit_chain}/1-{len_ref_align} DE [subseq from] mol:protein length:{len_ref_align}\n\n"
    # actual alignment
    formatted_alignment += f"{hit_id}_{hit_chain}/1-{len_ref_align}           "
    formatted_alignment += (
        "".join(
            [
                aa if query_align[i] != "-" else aa.lower()
                for i, aa in enumerate(ref_align)
            ]
        )
        + "\n"
    )
    formatted_alignment += "\n"

    return formatted_alignment


def fix_mmcif(path, chains, sequences, revision_date):
    """
    ensures that the mmcif header contains all
    necessary info for alphafold to correctly digest it
        * latest deposition date
        * which molecules are part of peptides
    """
    with open(path, "r") as infile:
        pdb_data = infile.readlines()

    # chains, amino acids, 3-letter codes, n
    # _entity_poly_seq.entity_id
    # _entity_poly_seq.num
    # _entity_poly_seq.mon_id
    # _entity_poly_seq.hetero
    # 1 1   GLN n
    # 1 2   SER n
    # ...
    # 2 1   THR n
    # ...
    for c, chain in reversed(list(enumerate(chains))):
        pdb_data.insert(
            2,
            "\n".join(
                [
                    f"{c+1 : <2} {r+1 : <4} {seq3(res).upper()} n"
                    if seq3(res).upper() != "XAA"
                    else f"{c+1 : <2} {r+1 : <4} UNK n"
                    for r, res in enumerate(sequences[c])
                ]
            )
            + "\n",
        )
    pdb_data.insert(
        2,
        "loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n_entity_poly_seq.hetero\n",
    )

    # chain info
    # loop_
    # _struct_asym.id
    # _struct_asym.entity_id
    # _struct_asym.pdbx_blank_PDB_chainid_flag
    # A 1 N
    # B 2 N
    pdb_data.insert(
        2,
        "\n".join([f"{chain} {c+1} N" for c, chain in enumerate(chains)])
        + "\n#\n",
    )
    pdb_data.insert(
        2,
        "loop_\n_struct_asym.id\n_struct_asym.entity_id\n_struct_asym.pdbx_blank_PDB_chainid_flag\n",
    )

    # peptide records
    pdb_data.insert(
        2,
        "\n".join(
            [
                f"{restype.upper()} peptide ? ? ? . "
                for restype in list(IUPACData.protein_letters_3to1.keys()) + ["Unk"]
            ]
        )
        + "\n#\n",
    )
    pdb_data.insert(
        2,
        "loop_\n_chem_comp.id\n_chem_comp.type\n_chem_comp.name\n_chem_comp.formula\n_chem_comp.formula_weight\n_chem_comp.mon_nstd_flag\n",
    )

    # the revision date is needed for AF to exclude templates newer than a given date
    pdb_data.insert(2, f"{revision_date}\n")
    pdb_data.insert(2, "_entry.id   pdb\n_pdbx_audit_revision_history.revision_date\n")

    asym_id_col = -1
    auth_id_col = -1
    atom_col_counter = 0
    for i, line in enumerate(pdb_data):
        if line.startswith("_atom_site"):
            if "label_asym_id" in line:
                asym_id_col = atom_col_counter
            elif "auth_asym_id" in line:
                auth_id_col = atom_col_counter
            atom_col_counter += 1

        if line.startswith("ATOM") and asym_id_col != auth_id_col:
            split_line = line.split(" ")
            column_indexes = [count for element, count in zip(split_line, range(len(split_line))) if element]
            split_line[column_indexes[asym_id_col]] = split_line[column_indexes[auth_id_col]]
            pdb_data[i] = " ".join(split_line)
    with open(path, "w") as out:
        out.write("".join(pdb_data))


def do_align(ref_seq, ref_model, query_seq, query_model, alignment_type="blast"):
    alignment = []
    if alignment_type == "blast":  # ref and query are sequences rather than structures
        aligner = Align.PairwiseAligner()
        aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -0.5
        #aligner.extend_gap_score = -0.5
        aln = aligner.align(ref_seq, query_seq)[0]
        try:  # compatibility between versions of Biopython
            alignment.append(aln[0])
            alignment.append(aln[1])
        except:
            alignment.append(aln.format().split("\n")[0])
            alignment.append(aln.format().split("\n")[2])
        print("\n".join(alignment))
    else:  # work with structural alignments instead
        ref_tempfile = tempfile.NamedTemporaryFile()
        query_tempfile = tempfile.NamedTemporaryFile()
        io = PDBIO()
        io.set_structure(ref_model)
        io.save(ref_tempfile.name)
        io.set_structure(query_model)
        io.save(query_tempfile.name)
        if alignment_type == "lddt_align":
            try:
                # save reference and query to disk on temporary file
                lddt_output = subprocess.check_output(
                    f"lDDT_align {ref_tempfile.name} {query_tempfile.name} --alignment_type horizontal",
                    shell=True,
                ).decode()
                alignment = re.findall(r"^[-A-Z]+$", lddt_output, flags=re.MULTILINE)
            except:
                print("Something went wrong when attempting to run lDDT_align")
                print(traceback.format_exc())
                sys.exit(1)
        elif alignment_type == "tmalign":
            try:
                # save reference and query to disk on temporary file
                tm_output = subprocess.check_output(
                    f"TMalign {ref_tempfile.name} {query_tempfile.name}", shell=True
                ).decode()
                alignment = re.findall(r"^[-A-Z]+$", tm_output, flags=re.MULTILINE)
            except:
                print("Something went wrong when attempting to run TMalign")
                print(traceback.format_exc())
                sys.exit(1)
        ref_tempfile.close()
        query_tempfile.close()
    return alignment


def get_target_data(paths, chains=None, is_fasta=False, is_mmcif=False):
    if not is_fasta:
        target_models = [load_PDB(target, is_mmcif=is_mmcif) for target in paths]
        if len(target_models) > 1:
            target_chains = (
                chains if chains else [c.id for model in target_models for c in model]
            )
            assert len(target_models) == len(
                target_chains
            ), "When superimposing chains from separate PDB files, these need to contain a single chain each, or you must specify one chain from each file with '--target_chains'"
            target_sequences = [
                get_fastaseq(model, chain)
                for model, chain in zip(target_models, target_chains)
            ]

        else:
            target_chains = chains if chains else [c.id for c in target_models[0]]
            target_sequences = [
                get_fastaseq(target_models[0], chain) for chain in target_chains
            ]
            # replicate the same model as many times as there are chains, so that all lists can be zipped later
            target_models = [pickle.loads(pickle.dumps(target_models[0], -1)) for i in range(len(target_chains))]
    else:  # fasta file containing one sequence per chain
        target_sequences = [record.seq for record in SeqIO.parse(paths[0], "fasta")]
        target_models = [None for s in target_sequences]
        target_chains = (
            chains
            if chains
            else [ascii_upperlower[i] for i, s in enumerate(target_sequences)]
        )

    return target_chains, target_sequences, target_models


def get_next_id(path):
    # these should go from 0000 to 9999, will complain if actual pdb ids are found e.g. 1abc
    cifs_in_path = sorted([Path(f).stem for f in Path(path).glob("*.cif")])
    if not cifs_in_path:
        return "0000"

    last_cif_id = cifs_in_path[-1]
    try:
        last_cif_id_int = int(last_cif_id)
    except:
        print(
            "ERROR: the template folder contains cif files not in the [0000-9999].cif format"
        )
        sys.exit(1)
    next_cif_id = last_cif_id_int + 1
    # returns e.g. 0004 if the latest file was 0003
    return str(next_cif_id).zfill(4)


def superimpose(ref_model, ref_chains, query_models, query_chains, alignment_type="tmalign"):
    backbone_atoms = ["CA", "C", "N", "O"]
    superimposer = Superimposer()
    for i, (ref_chain, query_chain) in enumerate(zip(ref_chains, query_chains)):
        query_atoms = []
        ref_atoms = []
        residues_to_delete = []
        ref_residues = ref_model[ref_chain].get_residues()

        # query models is a list of one model per chain every time, so we access the ith model and extract the ith chain
        query_model = query_models[i]
        query_residues = query_model[query_chain].get_residues()
        # align chains, get aligned atoms from reference and query
        alignment = do_align(
            "_",
            ref_model[ref_chain],
            "_",
            query_model[query_chain],
            alignment_type=alignment_type,
        )

        for ref_letter, query_letter in zip(alignment[0], alignment[1]):
            ref_aa = None
            query_aa = None
            if ref_letter != "-":
                ref_aa = next(ref_residues)
            if query_letter != "-":
                query_aa = next(query_residues)
            if ref_aa and query_aa:
                ref_ids = [atom.id for atom in ref_aa.get_atoms()]
                query_ids = [atom.id for atom in query_aa.get_atoms()]
                ref_atoms += [
                    atom for atom in ref_aa.get_atoms() if atom.id in backbone_atoms and atom.id in query_ids
                ]
                query_atoms += [
                    atom for atom in query_aa.get_atoms() if atom.id in backbone_atoms and atom.id in ref_ids
                ]
            elif query_aa:
                residues_to_delete.append(query_aa.get_full_id())
    
        #for _, _, chain, res in residues_to_delete:
        #    query_model[query_chain].detach_child(res)

        # superimpose queries, chain by chain
        superimposer.set_atoms(ref_atoms, query_atoms)
        superimposer.apply(query_model.get_atoms())

    if (
        query_models[0] is not query_models[1]
    ):  # we have superimposed chains from different models
        # merge query models by picking the right chains and adding them to the first
        for i, chain in enumerate(query_chains):
            current_chain = query_models[i][chain]
            for ch in query_models[i]:
                if ch is not current_chain:
                    query_models[i].detach_child(ch.id)
            current_chain.id = ascii_upperlower[i]
            if i > 0:
                if current_chain.id in query_models[0]:
                    query_models[0].detach_child(current_chain.id)
                query_models[0].add(current_chain)

    superimposed_model = query_models[0]

    return superimposed_model


def main():
    args = parse_args()
    fasta_target = is_fasta(args.target[0])

    if fasta_target:
        args.out_dir = f"{args.out_dir}/{Path(args.target[0]).stem}"

    # start generating output directory tree
    mmcif_path = Path(args.out_dir, "template_data", "mmcif_files")
    mmcif_path.mkdir(parents=True, exist_ok=True)
    # we always start from a  fake PDB id "0000", unless there are already templates from a previous run that we would like to add to
    next_id = get_next_id(mmcif_path) if args.append else "0000"

    args.align = args.align or args.align_tool

    if args.align and not args.align_tool:
        if fasta_target:
            args.align_tool = "blast"
        else:
            args.align_tool = "tmalign"

    # load target data if needed
    if args.align or args.superimpose:
        target_chains, target_sequences, target_models = get_target_data(
            args.target,
            chains=args.target_chains,
            is_fasta=fasta_target,
            is_mmcif=args.mmcif_target,
        )

    # Handling the template file: convert to a compatible mmCIF file, write sequences to pdb_seqres.txt
    template_model = load_PDB(args.template, is_mmcif=args.mmcif_template)
    # template sequences are needed to write out the pdb_seqres.txt file
    template_chains = (
        args.template_chains if args.template_chains else [c.id for c in template_model]
    )
    remove_extra_chains(template_model, template_chains)
    remove_hetatms(template_model)
    template_sequences = [
        get_fastaseq(template_model, chain) for chain in template_chains
    ]

    io = MMCIFIO()
    template_mmcif_path = os.path.join(
        args.out_dir, "template_data", "mmcif_files", f"{next_id}.cif"
    )

    if args.superimpose:  # modify template
        # superimpose target chains to template, then save those as template mmcif, and realign to itself
        target_model = superimpose(
            template_model, template_chains, target_models, target_chains, alignment_type=args.align_tool
        )
        template_model = target_model
        template_sequences = target_sequences
        template_chains = target_chains

    if args.inpaint_clashes:
        print("Deleting clashing residues...", end=" ")
        n_deleted, template_model = detect_and_remove_clashes(template_model)
        template_sequences = [
            get_fastaseq(template_model, chain) for chain in template_chains
        ]
        print(f"{n_deleted} clashes found.")

    io.set_structure(template_model)
    io.save(template_mmcif_path)

    fix_mmcif(
        template_mmcif_path, template_chains, template_sequences, args.revision_date
    )

    pdb_seqres_path = Path(args.out_dir, "template_data", "pdb_seqres.txt").resolve()
    write_seqres(
        pdb_seqres_path,
        template_sequences,
        template_chains,
        seq_id=next_id,
        append=args.append,
    )

    # extra flagfile for AF usage
    af_flagfile_path = Path(args.out_dir, "template_data", "templates.flag")
    if not af_flagfile_path.is_file():  # don't overwrite file if already there
        with open(af_flagfile_path, "w") as flagfile:
            flagfile.write(f"--template_mmcif_dir={mmcif_path.resolve()}\n")
            flagfile.write(f"--pdb_seqres_database_path={pdb_seqres_path}\n")
            if args.align:  # means we are not going to let AF overwrite pdb_hits.sto
                flagfile.write("--use_precomputed_msas\n")
    """
    Handling targets: here alignments are performed against the template
    The target can either be one or more PDB files, or a fasta file containing sequences.
    
    If a fasta file is submitted, then sequence alignments will be performed
    If one or more PDBs are submitted, then either sequence or structural alignments can be performed
    
    If multiple model PDBs are submitted, then we are superimposing several unbound chains to the same template
    """
    if (
        args.align
    ):  # only if an alignment tool is selected, otherwise leave it to AlphaFold's template search

        assert len(target_chains) == len(
            template_chains
        ), f"The number of chains to align from target ({target_chains}) doesn't match the number of chains in the template ({template_chains}). Make sure that the files contain the same number of chains or select the chains that should be paired with --target_chains, --template_chains"
        for (
            i,
            (
                template_chain,
                template_sequence,
                target_chain,
                target_sequence,
                target_model,
            ),
        ) in enumerate(
            zip(
                template_chains,
                template_sequences,
                target_chains,
                target_sequences,
                target_models,
            )
        ):
            msa_chain = ascii_upperlower[i]
            this_template_model = pickle.loads(pickle.dumps(template_model, -1))
            this_target_model = pickle.loads(pickle.dumps(target_model, -1))
            if not fasta_target:
                remove_extra_chains(this_template_model, [template_chain])
                remove_extra_chains(this_target_model, [target_chain])
            print(f"\nAligning fasta sequence {i+1} (seq: {target_sequence[0:10]}...) to template chain {template_chain} (seq: {template_sequence[0:10]}...)")
            alignment = do_align(
                template_sequence,
                this_template_model,
                target_sequence,
                this_target_model,
                alignment_type=args.align_tool,
            )
            sto_alignment = format_alignment_stockholm(
                alignment, hit_id=next_id, hit_chain=template_chain
            )
            
            
            msa_path = f"msas/{msa_chain}"
            
            # write alignment to file
            Path(args.out_dir, msa_path).mkdir(parents=True, exist_ok=True)
            with open(
                Path(args.out_dir, msa_path, "pdb_hits.sto"),
                mode="a" if args.append else "w",
            ) as pdb_hits:
                for line in sto_alignment:
                    pdb_hits.write(line)

    if not fasta_target:
        print(
            f"Run AlphaFold with, e.g.:\npython run_alphafold.py --fasta_paths target.fasta --flagfile databases.flag --flagfile {af_flagfile_path} --output_dir {Path(args.out_dir).parents[0]} --cross_chain_templates --dropout --model_preset='multimer_v2' --separate_homomer_msas"
        )
        print(
            "*** NB: the name of the fasta target should be the same as the name of the folder containing the output msas: (e.g.  if the fasta target file is 'target.fasta', then --output_dir='somedir/target' ***"
        )
    else:
        print(
            f"Run AlphaFold with, e.g.:\npython run_alphafold.py --fasta_paths {args.target[0]} --flagfile databases.flag --flagfile {af_flagfile_path} --output_dir {Path(args.out_dir).parents[0]} --cross_chain_templates --dropout --model_preset='multimer_v2' --separate_homomer_msas"
        )


if __name__ == "__main__":
    main()
