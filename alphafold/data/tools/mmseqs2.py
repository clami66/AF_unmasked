import math
import shutil
import subprocess
from pathlib import Path
from absl import logging
from typing import List, Union
from alphafold.data.tools import utils

logging.set_verbosity(logging.INFO)
logging.use_absl_handler()


class MMseqs2:
    def __init__(self,
                 uniref_db: str,
                 binary_path: str,
                 metagenomic_db: str = None,
                 expand_eval: float = math.inf,
                 diff: int = 3000,
                 s: float = 8,
                 db_load_mode: int = 0,
                 n_cpu: int = 32,
                 gpu: bool = True,
                 gpu_server: bool = False,
                 msa_out_dir: Path = None,
                 max_accept: int = 100000,
                 ):
        self.mmseqs = Path(binary_path)
        self.uniref_db = Path(uniref_db)
        self.metagenomic_db = Path(metagenomic_db) if metagenomic_db is not None else None
        self.expand_eval = str(expand_eval)
        self.diff = str(diff)
        self.s = s
        self.db_load_mode = str(db_load_mode)
        self.threads = str(n_cpu)
        self.gpu = gpu
        self.gpu_server = gpu_server
        self.base = msa_out_dir

        self.align_eval = "10"
        self.qsc = "0.8"
        self.max_accept = str(max_accept)

        if not self.uniref_db.with_suffix(".dbtype").is_file():
            raise FileNotFoundError(f"Database {uniref_db} does not exist")
        if (
            not self.uniref_db.with_suffix(".idx").is_file()
            and not self.uniref_db.with_suffix(".idx.index").is_file()
        ):
            logging.info("Uniref search does not use index")
            self.db_load_mode = "0"
            dbSuffix1 = "_seq"
            dbSuffix2 = "_aln"
        else:
            logging.info("Uniref search with index")
            dbSuffix1 = ".idx"
            dbSuffix2 = ".idx"

        self.uniref_db_1 = Path(uniref_db + dbSuffix1)
        self.uniref_db_2 = Path(uniref_db + dbSuffix2)

        if self.metagenomic_db is not None:
            if not self.metagenomic_db.with_suffix(".dbtype").is_file():
                raise FileNotFoundError(f"Database {metagenomic_db} does not exist")
            if (
                not self.metagenomic_db.with_suffix(".idx").is_file()
                and not self.metagenomic_db.with_suffix(".idx.index").is_file()
            ):
                logging.info("Metagenomic DB search does not use index")
                dbSuffix1 = "_seq"
                dbSuffix2 = "_aln"
            else:
                logging.info("Metagenomic DB search with index")
                dbSuffix1 = ".idx"
                dbSuffix2 = ".idx"
            self.metagenomic_db_1 = Path(str(self.metagenomic_db) + dbSuffix1)
            self.metagenomic_db_2 = Path(str(self.metagenomic_db) + dbSuffix2)

        self.search_param = ["--num-iterations", "3",
                             "--db-load-mode", self.db_load_mode,
                             "-a",
                             "-e", "0.1",
                             "--max-seqs", "10000"]

        if self.gpu:
            # gpu version only supports ungapped prefilter currently
            self.search_param.extend(["--gpu", "1", "--prefilter-mode", "1"])
        else:
            # self.search_param.extend(["--prefilter-mode", "0"])
             # sensitivy can only be set for non-gpu version, gpu version runs at max sensitivity
            if self.s is not None:
                self.search_param += ["-s", "{:.1f}".format(self.s)]
            else:
                self.search_param += ["--k-score", "'seq:96,prof:80'"]
        if gpu_server:
            self.search_param += ["--gpu-server", "1"]

        self.filter_param = ["--filter-msa", "1",
                             "--filter-min-enable", "1000",
                             "--diff", self.diff,
                             "--qid", "0.0,0.2,0.4,0.6,0.8,1.0",
                             "--qsc", "0",
                             "--max-seq-id", "0.95"]
        self.expand_param = ["--expansion-mode", "0",
                             "-e", self.expand_eval,
                             "--expand-filter-clusters", "1",
                             "--max-seq-id", "0.95"]

    def run_mmseqs(self, params: List[Union[str, Path]]):
        params_log = " ".join(str(i) for i in params)
        logging.info(f"Running {self.mmseqs} {params_log}")
        subprocess.check_call([self.mmseqs] + params)


    def run_mmseqs_parallel(self, params_sets):
        p1 = subprocess.Popen([self.mmseqs] + params_sets[0]) 
        p3 = subprocess.Popen([self.mmseqs] + params_sets[2])
        p1.wait()
        p2 = subprocess.Popen([self.mmseqs] + params_sets[1])
        p2.wait()
        p3.wait()
        """
        params_log_sets = [" ".join(str(i) for i in params) for params in params_sets]
        for params_log in params_log_sets:
            logging.info(f"Running in parallel {self.mmseqs} {params_log}")
        for p in procs:
            p.wait()
        """

    def cleanup(self, base):
        self.run_mmseqs(["rmdb", base.joinpath("res_exp_realign_filter")])
        self.run_mmseqs(["rmdb", base.joinpath("res_exp_realign")])
        self.run_mmseqs(["rmdb", base.joinpath("res_exp")])
        self.run_mmseqs(["rmdb", base.joinpath("res")])
        self.run_mmseqs(["rmdb", base.joinpath("uniref.a3m")])
        self.run_mmseqs(["rmdb", base.joinpath("final.a3m")])
        self.run_mmseqs(["rmdb", base.joinpath("prof_res")])
        self.run_mmseqs(["rmdb", base.joinpath("prof_res_h")])
        self.run_mmseqs(["rmdb", base.joinpath("qdb")])
        self.run_mmseqs(["rmdb", base.joinpath("qdb_h")])
        shutil.rmtree(base.joinpath("tmp"))

        if self.metagenomic_db is not None:
            self.run_mmseqs(["rmdb", base.joinpath("res_env_exp_realign_filter")])
            self.run_mmseqs(["rmdb", base.joinpath("res_env_exp_realign")])
            self.run_mmseqs(["rmdb", base.joinpath("res_env_exp")])
            self.run_mmseqs(["rmdb", base.joinpath("res_env")])
            self.run_mmseqs(["rmdb", base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m")])
            shutil.rmtree(base.joinpath("tmp_env"))

    def query(
        self, fasta_path: str, max_sequences:int = None
    ):

        delete = False if self.base is not None else True
        with utils.tmpdir_manager(base_dir=self.base, prefix="alignments", delete=delete) as query_tmp_dir:
            base = Path(query_tmp_dir)
            self.run_mmseqs(
                ["createdb", fasta_path,  base.joinpath("qdb"), "--shuffle", "0"],
            )

            self.run_mmseqs(["search",
                             base.joinpath("qdb"),
                             self.uniref_db,
                             base.joinpath("res"),
                             base.joinpath("tmp"),
                             "--threads", self.threads] + self.search_param)
            self.run_mmseqs(["mvdb", base.joinpath("tmp/latest/profile_1"), base.joinpath("prof_res")])
            self.run_mmseqs(["lndb", base.joinpath("qdb_h"), base.joinpath("prof_res_h")])

            if self.metagenomic_db is None:
                self.run_mmseqs(["expandaln",
                                 base.joinpath("qdb"),
                                 self.uniref_db_1,
                                 base.joinpath("res"),
                                 self.uniref_db_2,
                                 base.joinpath("res_exp"),
                                 "--db-load-mode", self.db_load_mode,
                                 "--threads", self.threads] + self.expand_param)
                self.run_mmseqs(["align",
                                base.joinpath("prof_res"),
                                self.uniref_db_1,
                                base.joinpath("res_exp"),
                                base.joinpath("res_exp_realign"),
                                "--db-load-mode", self.db_load_mode,
                                "-e", self.align_eval,
                                "--max-accept", self.max_accept,
                                "--threads", self.threads,
                                "--alt-ali", "10", "-a"])
            else:
                self.run_mmseqs_parallel([
                                ["expandaln",
                                 base.joinpath("qdb"),
                                 self.uniref_db_1,
                                 base.joinpath("res"),
                                 self.uniref_db_2,
                                 base.joinpath("res_exp"),
                                 "--db-load-mode", self.db_load_mode,
                                 "--threads", str(int(self.threads)//2)] + self.expand_param,
                                ["align",
                                base.joinpath("prof_res"),
                                self.uniref_db_1,
                                base.joinpath("res_exp"),
                                base.joinpath("res_exp_realign"),
                                "--db-load-mode", self.db_load_mode,
                                "-e", self.align_eval,
                                "--max-accept", self.max_accept,
                                "--threads", str(int(self.threads)//2),
                                "--alt-ali", "10", "-a"],
                                ["search",
                                 base.joinpath("prof_res"),
                                 self.metagenomic_db,
                                 base.joinpath("res_env"),
                                 base.joinpath("tmp_env"),
                                 "--threads", str(int(self.threads)//2)] + self.search_param ])

            self.run_mmseqs(["filterresult", base.joinpath("qdb"),
                             self.uniref_db_1,
                             base.joinpath("res_exp_realign"),
                             base.joinpath("res_exp_realign_filter"),
                             "--db-load-mode", self.db_load_mode,
                             "--qid", "0",
                             "--qsc", self.qsc,
                             "--diff", "0",
                             "--threads", self.threads,
                             "--max-seq-id", "1.0",
                             "--filter-min-enable", "100"])
            self.run_mmseqs(["result2msa",
                             base.joinpath("qdb"),
                             self.uniref_db_1,
                             base.joinpath("res_exp_realign_filter"),
                             base.joinpath("uniref.a3m"),
                             "--msa-format-mode", "3",
                             "--db-load-mode", self.db_load_mode,
                             "--threads", self.threads] + self.filter_param)

            if self.metagenomic_db is not None:
                """
                self.run_mmseqs(["search",
                                 base.joinpath("prof_res"),
                                 self.metagenomic_db,
                                 base.joinpath("res_env"),
                                 base.joinpath("tmp_env"),
                                 "--threads", self.threads] + self.search_param)
                """
                self.run_mmseqs(["expandaln",
                                 base.joinpath("prof_res"),
                                 self.metagenomic_db_1,
                                 base.joinpath("res_env"),
                                 self.metagenomic_db_2,
                                 base.joinpath("res_env_exp"),
                                 "-e", self.expand_eval,
                                 "--expansion-mode", "0",
                                 "--db-load-mode", self.db_load_mode,
                                 "--threads", self.threads])
                self.run_mmseqs(["align",
                                 base.joinpath("tmp_env/latest/profile_1"),
                                 self.metagenomic_db_1,
                                 base.joinpath("res_env_exp"),
                                 base.joinpath("res_env_exp_realign"),
                                 "--db-load-mode", self.db_load_mode,
                                 "-e", self.align_eval,
                                 "--max-accept", self.max_accept, ############# test this
                                 "--threads", self.threads,
                                 "--alt-ali", "10", "-a"])
                self.run_mmseqs(["filterresult",
                                 base.joinpath("qdb"),
                                 self.metagenomic_db_1,
                                 base.joinpath("res_env_exp_realign"),
                                 base.joinpath("res_env_exp_realign_filter"),
                                 "--db-load-mode", self.db_load_mode,
                                 "--qid", "0",
                                 "--qsc", self.qsc,
                                 "--diff", "0",
                                 "--max-seq-id", "1.0",
                                 "--threads", self.threads,
                                 "--filter-min-enable", "100"])
                self.run_mmseqs(["result2msa",
                                 base.joinpath("qdb"),
                                 self.metagenomic_db_1,
                                 base.joinpath("res_env_exp_realign_filter"),
                                 base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m"),
                                 "--msa-format-mode", "3",
                                 "--db-load-mode", self.db_load_mode,
                                 "--threads", self.threads] + self.filter_param)

                self.run_mmseqs(["mergedbs", base.joinpath("qdb"), base.joinpath("final.a3m"), base.joinpath("uniref.a3m"), base.joinpath("bfd.mgnify30.metaeuk30.smag30.a3m")])
            else:
                self.run_mmseqs(["mvdb", base.joinpath("uniref.a3m"), base.joinpath("final.a3m")])

            self.run_mmseqs(["unpackdb", base.joinpath("final.a3m"), base.joinpath("."), "--unpack-name-mode", "0", "--unpack-suffix", ".a3m"])

            self.cleanup(base)

            with open(base.joinpath("0.a3m")) as f:
                a3m = f.readlines()
            a3m = "".join([line for line in a3m if line.strip() and not line.startswith("#")])
            
            #if max_sequences is not None:
            #    a3m = "\n".join(a3m.split("\n")[:max_sequences*2])
        raw_output = dict(
            a3m=a3m,
            output=None,
            stderr=None,
            n_iter=3,
            s=float(self.s),
            e_value=float(self.align_eval))
        return [raw_output]
