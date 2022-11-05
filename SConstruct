#!/usr/bin/env python
# -*- coding: utf-8 -*-
import nestly
from nestly import scons as nestly_scons
import os
from SCons.Script import Environment
import SCons.Script as Script
import sys

#### Set up command line arguments/options

all_actions = ['run_partis', 'run_linearham', 'build_partis_linearham']

# partis arguments

Script.AddOption("--run-partis",
        dest="run_partis",
        action="store_true",
        default=False,
        help="Should we run partis?")

Script.AddOption("--fasta-path",
        dest="fasta_path",
        type="str",
        default=None,
        help="The repertoire FASTA file path.")

Script.AddOption("--all-clonal-seqs",
        dest="all_clonal_seqs",
        action="store_true",
        default=False,
        help="Should we assume all the sequences in the FASTA file are clonal?")

Script.AddOption("--locus",
        dest="locus",
        type="str",
        default="igh",
        help="Which immunoglobulin locus are we doing inference on?")

# linearham arguments

Script.AddOption("--run-linearham",
        dest="run_linearham",
        action="store_true",
        default=False,
        help="Should we run linearham?")

Script.AddOption("--partition-index",
        dest="partition_index",
        type="str",
        default=None,
        help="Zero-based index of partition from which to choose a cluster in the partis list of most likely partitions (default: index of most likely partition). Anything passed will be interpreted as the index of the partition step to be used as in scripts/parse_cluter.py.")

Script.AddOption("--cluster-index",
        dest="cluster_index",
        type="str",
        default=None,
        help="Zero-based index of cluster to use in the partition specified by --partition-index (default: 0). Anything passed will be interpreted as the index of the cluster to be used as in scripts/parse_cluter.py.")

Script.AddOption("--cluster-seed-unique-id",
        dest="cluster_seed_unique_id",
        type="str",
        default=None,
        help="A string specifying unique id of the partis seed sequence (see partis --seed-unique-id https://github.com/psathyrella/partis/blob/master/docs/subcommands.md#--seed-unique-id-id). The uid passed will be used to parse the cluster(s) containing this sequence in the partition specified by --partition-index as in scripts/parse_cluter.py.")

Script.AddOption("--template-path",
        dest="template_path",
        type="str",
        default="templates/revbayes_template.rev",
        help="The Rev template path.")

Script.AddOption("--mcmc-iter",
        dest="mcmc_iter",
        type="str",
        default="10000",
        help="How many RevBayes MCMC iterations should we use?")

Script.AddOption("--mcmc-thin",
        dest="mcmc_thin",
        type="str",
        default="10",
        help="What RevBayes MCMC thinning frequency should we use?")

Script.AddOption("--tune-iter",
        dest="tune_iter",
        type="str",
        default="5000",
        help="How many RevBayes tuning iterations should we use?")

Script.AddOption("--tune-thin",
        dest="tune_thin",
        type="str",
        default="100",
        help="What RevBayes tuning thinning frequency should we use?")

Script.AddOption("--num-rates",
        dest="num_rates",
        type="str",
        default="4",
        help="The number of gamma rate categories.")

Script.AddOption("--burnin-frac",
        dest="burnin_frac",
        type="str",
        default="0.1",
        help="What fraction of MCMC burnin should we use?")

Script.AddOption("--subsamp-frac",
        dest="subsamp_frac",
        type="str",
        default="0.05",
        help="What bootstrap sampling fraction should we use?")

Script.AddOption("--rng-seed",
        dest="rng_seed",
        type="str",
        default="0",
        help="The RNG seed.")

Script.AddOption("--lineage-unique-ids",
        dest="lineage_unique_ids",
        type="str",
        default=None,
        help=",-separated list of names of the sequence(s) of interest for analyzing lineages. ")

Script.AddOption("--asr-pfilters",
        dest="asr_pfilters",
        type="str",
        default="0.1",
        help="The ancestral sequence posterior probability threshold.")

Script.AddOption("--partis-yaml-file",
        dest="partis_yaml_file",
        type="str",
        default=None,
        help="An optional partis output YAML file.")

# partis/linearham arguments

Script.AddOption("--build-partis-linearham",
        dest="build_partis_linearham",
        action="store_true",
        default=False,
        help="Should we build partis and linearham?")

Script.AddOption("--parameter-dir",
        dest="parameter_dir",
        type="str",
        default=None,
        help="An optional directory of partis parameter files.")

Script.AddOption("--outdir",
        dest="outdir",
        default="output",
        help="The output directory.")


#### Process command line arguments/options

def process_multiarg(arg_str, type, delim):
    return [type(arg) for arg in arg_str.split(delim)]

def get_options(env):
    return dict(
        # partis arguments
        run_partis = env.GetOption("run_partis"),
        fasta_path = env.GetOption("fasta_path"),
        all_clonal_seqs = env.GetOption("all_clonal_seqs"),
        locus = env.GetOption("locus"),

        # linearham arguments
        run_linearham = env.GetOption("run_linearham"),
        cluster_seed_unique_id = str(env.GetOption("cluster_seed_unique_id")) if env.GetOption("cluster_seed_unique_id") is not None else None,
        cluster_index = int(env.GetOption("cluster_index")) if env.GetOption("cluster_index") is not None else None,
        partition_index = int(env.GetOption("partition_index")) if env.GetOption("partition_index") is not None else None,
        template_path = env.GetOption("template_path"),
        mcmc_iter = process_multiarg(env.GetOption("mcmc_iter"), int, ","),
        mcmc_thin = process_multiarg(env.GetOption("mcmc_thin"), int, ","),
        tune_iter = process_multiarg(env.GetOption("tune_iter"), int, ","),
        tune_thin = process_multiarg(env.GetOption("tune_thin"), int, ","),
        num_rates = process_multiarg(env.GetOption("num_rates"), int, ","),
        burnin_frac = process_multiarg(env.GetOption("burnin_frac"), float, ","),
        subsamp_frac = process_multiarg(env.GetOption("subsamp_frac"), float, ","),
        rng_seed = process_multiarg(env.GetOption("rng_seed"), int, ","),
        lineage_unique_ids = process_multiarg(env.GetOption("lineage_unique_ids"), str, ",") if env.GetOption("lineage_unique_ids") is not None else None,
        asr_pfilters = process_multiarg(env.GetOption("asr_pfilters"), float, ","),
        partis_yaml_file = env.GetOption("partis_yaml_file"),

        # partis/linearham arguments
        build_partis_linearham = env.GetOption("build_partis_linearham"),

        # general arguments
        parameter_dir = env.GetOption("parameter_dir"),
        outdir = env.GetOption("outdir")
    )

env = Environment(ENV = os.environ)
options = get_options(env)
if not options["build_partis_linearham"]:
    env.SConsignFile(os.path.join(options["outdir"], ".sconsign"))

if not any(options[a] for a in all_actions) and '--help' not in sys.argv:
    raise Exception('No action specified. Choose from %s' % ', '.join('--' + a.replace('_', '-') for a in all_actions))

lhdir = os.getcwd()  # os.path.dirname(os.path.realpath(__file__))  # would like to do it this way, but __file__ isn't set for scripts, only modules
if os.path.abspath(options["outdir"]).find(lhdir) != 0:
    raise Exception("in order to avoid scons getting confused and not doing anything, --outdir has to be a subdir of the main linearham dir (%s), but it was set to %s" % (lhdir, os.path.abspath(options["outdir"])))

#### Set up the nesting structure

nest = nestly_scons.SConsWrap(nestly.Nest(), dest_dir=options["outdir"], alias_environment=env)

def default_label(d):
    return d.get("id")

#### Record git commit and tag

@nest.add_target()
def commit_and_tag(outdir, c):
    log = env.Command(os.path.join(outdir, "git.log"), "", "git rev-parse --verify HEAD >> $TARGET && git describe --dirty --always >> $TARGET")
    env.AlwaysBuild(log)
    return log

#### Install partis and linearham (if necessary)

if options["build_partis_linearham"]:

    @nest.add_target()
    def partis_build(outdir, c):
        partisbuild = env.Command("lib/partis/packages/ham/bcrham", "",
                           "cd lib/partis && ./bin/build.sh")
        env.AlwaysBuild(partisbuild)
        return partisbuild

    @nest.add_target()
    def linearham_build(outdir, c):
        libptpll_env = env.Clone()
        libptpll_env.VariantDir("_build/libptpll", "lib/libptpll", duplicate=0)
        libptpll_env.Command("_build/libptpll/libptpll_static.a", "",
                             "cp -r lib/libptpll _build/ && cd _build/libptpll && " + \
                             "make && cp -t . _build/src/libptpll_static.a " + \
                             "_build/lib/_prefix/lib/libpll_algorithm.a " + \
                             "_build/lib/_prefix/lib/libpll_optimize.a " + \
                             "_build/lib/_prefix/lib/libpll_tree.a " + \
                             "_build/lib/_prefix/lib/libpll_util.a " + \
                             "_build/lib/_prefix/lib/libpll.a " + \
                             "_build/lib/lesplace/src/liblesplace-static.a")

        linearham_env = env.Clone()
        linearham_env.Append(CPPPATH=["lib/eigen", "lib/fast-cpp-csv-parser",
                                      "lib/libptpll/src", "lib/tclap/include",
                                      "_build/libptpll/_build/lib/_prefix/include"])
        linearham_env.Append(CCFLAGS=["-pthread", "-std=c++11", "-g"])
        linearham_env.Append(LIBPATH=["_build/libptpll"])
        linearham_env.Append(LIBS=["ptpll_static", "pll_algorithm", "pll_optimize",
                                   "pll_tree", "pll_util", "pll", "lesplace-static",
                                   "gsl", "blas", "yaml-cpp", "pthread"])
        linearham_env.Append(LINKFLAGS=["-g"])
        linearham_env.VariantDir("_build/linearham", "src")
        linearham_env.StaticLibrary(
            target="_build/linearham/liblinearham.a",
            source=Glob("_build/linearham/*.cpp",
                        exclude=["_build/linearham/linearham.cpp"])
        )
        linearham_env.Append(CPPPATH=["src"])
        linearham_env.Append(LIBPATH=["_build/linearham"])
        linearham_env.Prepend(LIBS=["linearham"])
        linearham_bin = linearham_env.Program(target="_build/linearham/linearham",
                                              source="_build/linearham/linearham.cpp")

        test_env = linearham_env.Clone()
        test_env.VariantDir("_build/test", "test")
        test_env.Program(target="_build/test/test", source="_build/test/test.cpp")
        return linearham_bin


#### Run partis (if necessary)

if options["run_partis"]:
    if options["fasta_path"] is None:
        raise Exception("--fasta-path is required for --run-partis")

    @nest.add_target()
    def partis_output(outdir, c):
        partis_mode = "annotate --all-seqs-simultaneous" if options["all_clonal_seqs"] else "partition"
# ----------------------------------------------------------------------------------------
        partis_parameter_dir = options["parameter_dir"].rstrip("/") + " --refuse-to-cache-parameters" if options["parameter_dir"] is not None else os.path.join(outdir, "parameter_dir")
        partis_output = env.Command(
            [os.path.join(outdir, filename) for filename in
                ["partis_run.yaml", "partis_run.stdout.log"]],
            options["fasta_path"],
            "lib/partis/bin/partis " + partis_mode \
                + " --infname $SOURCE" \
                + " --parameter-dir " + partis_parameter_dir \
                + " --locus " + options["locus"] \
                + " --extra-annotation-columns linearham-info" \
                + " --outfname ${TARGETS[0]} > ${TARGETS[1]}")
        env.Depends(partis_output, "lib/partis/packages/ham/bcrham")
        return partis_output


#### Run linearham (if necessary)

if options["run_linearham"]:

    @nest.add_target()
    def partis_yaml_file(outdir, c):
        default_outfname = os.path.join(outdir, "partis_run.yaml")
        if options["partis_yaml_file"] is not None:
# ----------------------------------------------------------------------------------------
            assert options["parameter_dir"] is not None, "Specify both --partis-yaml-file and --parameter-dir."
            if options["partis_yaml_file"] == default_outfname:
                # scons complains about dependency cycles if source and target have the same name
                default_outfname = os.path.join(outdir, "partis_run.lh.yaml")
            partis_yaml_file = env.Command(
                default_outfname,
                options["partis_yaml_file"],
                "lib/partis/bin/partis get-linearham-info" \
                    + " --outfname $SOURCE" \
                    + " --parameter-dir " + options["parameter_dir"] \
                    + " --linearham-info-fname $TARGET")
            env.Depends(partis_yaml_file, "lib/partis/packages/ham/bcrham")
            return str(partis_yaml_file[0])
        else:
            return default_outfname

    @nest.add_target()
    def hmm_param_dir(outdir, c):
        linearham_parameter_dir = options["parameter_dir"].rstrip("/") if options["parameter_dir"] is not None else os.path.join(outdir, "parameter_dir")
        return linearham_parameter_dir + "/hmm/hmms"

    @nest.add_nest(label_func=default_label)
    def cluster(c):
        cluster_info = {"cluster_index": options["cluster_index"], "partition_index": options["partition_index"], "cluster_seed_unique_id": options["cluster_seed_unique_id"]}
        if options["cluster_seed_unique_id"] is not None:
            cluster_info["id"] = "cluster-" + str(options["cluster_seed_unique_id"])
        elif options["cluster_index"] is not None:
            cluster_info["id"] = "cluster-" + str(options["cluster_index"])
        else:
            cluster_info["id"] = "cluster-0"
        return [cluster_info]

    if options["lineage_unique_ids"] is not None:
        @nest.add_nest(label_func=default_label)
        def lineage(c):
            return [{"id": "lineage_" + lid, "lineage_seq_unique_id": lid}
                    for lid in options["lineage_unique_ids"]]

    @nest.add_target()
    def _parse_cluster(outdir, c):
        cluster_fasta_file, cluster_yaml_file = env.Command(
            [os.path.join(outdir, outfile) for outfile in ["cluster_seqs.fasta", "cluster.yaml"]],
            c["partis_yaml_file"],
            "scripts/parse_cluster.py $SOURCE" \
                + " --indel-reversed-seqs" \
                + ((" --seed-unique-id " + str(c["cluster"]["cluster_seed_unique_id"])) if c["cluster"]["cluster_seed_unique_id"] is not None else "") \
                + ((" --cluster-index " + str(c["cluster"]["cluster_index"])) if c["cluster"]["cluster_index"] is not None else "") \
                + ((" --partition-index " + str(c["cluster"]["partition_index"])) if c["cluster"]["partition_index"] is not None else "") \
                + ((" --seed-seq " + c["lineage"]["lineage_seq_unique_id"]) if options["lineage_unique_ids"] is not None else "") \
# ----------------------------------------------------------------------------------------
                + ((" --glfo-dir " + os.path.join(options["parameter_dir"], "/hmm/germline-sets")) if os.path.splitext(c["partis_yaml_file"])[1] == '.csv' else "") \
                + ((" --locus " + options["locus"]) if os.path.splitext(c["partis_yaml_file"])[1] == '.csv' else "") \
                + " --fasta-output-file ${TARGETS[0]}" \
                + " --yaml-output-file ${TARGETS[1]}")
        env.Depends(cluster_fasta_file, "scripts/parse_cluster.py")
        return cluster_fasta_file, cluster_yaml_file

    @nest.add_target()
    def cluster_fasta_file(outdir, c):
        return c['_parse_cluster'][0]

    @nest.add_target()
    def cluster_yaml_file(outdir, c):
        return c['_parse_cluster'][1]

    @nest.add_nest(label_func=default_label)
    def revbayes_setting(c):
        return [{"id": "mcmciter" + str(mcmc_iter) + "_mcmcthin" + str(mcmc_thin) + \
                       "_tuneiter" + str(tune_iter) + "_tunethin" + str(tune_thin) + \
                       "_numrates" + str(num_rates) + "_rngseed" + str(rng_seed),
                 "mcmc_iter": mcmc_iter, "mcmc_thin": mcmc_thin,
                 "tune_iter": tune_iter, "tune_thin": tune_thin,
                 "num_rates": num_rates, "rng_seed": rng_seed}
                for mcmc_iter in options["mcmc_iter"]
                for mcmc_thin in options["mcmc_thin"]
                for tune_iter in options["tune_iter"]
                for tune_thin in options["tune_thin"]
                for num_rates in options["num_rates"]
                for rng_seed in options["rng_seed"]]

    @nest.add_target()
    def revbayes_rev_file(outdir, c):
        revbayes_rev_file = env.Command(
            os.path.join(outdir, "revbayes_run.rev"),
            [options["template_path"], c["cluster_fasta_file"]],
            "scripts/generate_revbayes_rev_file.py ${SOURCES[0]}" \
                + " --fasta-path ${SOURCES[1]}" \
                + " --mcmc-iter " + str(c["revbayes_setting"]["mcmc_iter"]) \
                + " --mcmc-thin " + str(c["revbayes_setting"]["mcmc_thin"]) \
                + " --tune-iter " + str(c["revbayes_setting"]["tune_iter"]) \
                + " --tune-thin " + str(c["revbayes_setting"]["tune_thin"]) \
                + " --num-rates " + str(c["revbayes_setting"]["num_rates"]) \
                + " --seed " + str(c["revbayes_setting"]["rng_seed"]) \
                + " --output-path $TARGET")
        env.Depends(revbayes_rev_file, "scripts/generate_revbayes_rev_file.py")
        return revbayes_rev_file

    @nest.add_target()
    def revbayes_output(outdir, c):
        revbayes_output = env.Command(
            [os.path.join(outdir, filename) for filename in
                ["revbayes_run.trees", "revbayes_run.log", "revbayes_run.stdout.log"]],
            c["revbayes_rev_file"],
            "lib/revbayes/projects/cmake/rb $SOURCE > ${TARGETS[2]}")
        env.Depends(revbayes_output, "lib/revbayes/projects/cmake/rb")
        return revbayes_output

    @nest.add_target()
    def linearham_intermediate_output(outdir, c):
        linearham_intermediate_output = env.Command(
            os.path.join(outdir, "lh_revbayes_run.trees"),
            [c["revbayes_output"][0], c["cluster_yaml_file"]],
            "_build/linearham/linearham --pipeline" \
                + " --yaml-path ${SOURCES[1]}" \
                # always use 0 here since we have created a partis yaml with only one cluster
                + " --cluster-ind 0" \
                + " --hmm-param-dir " + c["hmm_param_dir"] \
                + " --seed " + str(c["revbayes_setting"]["rng_seed"]) \
                + " --num-rates " + str(c["revbayes_setting"]["num_rates"]) \
                + " --input-path ${SOURCES[0]}" \
                + " --output-path $TARGET")
        env.Depends(linearham_intermediate_output, "_build/linearham/linearham")
        return linearham_intermediate_output

    @nest.add_nest(label_func=default_label)
    def linearham_setting(c):
        return [{"id": "burninfrac" + str(burnin_frac) + "_subsampfrac" + str(subsamp_frac),
                 "burnin_frac": burnin_frac, "subsamp_frac": subsamp_frac}
                for burnin_frac in options["burnin_frac"]
                for subsamp_frac in options["subsamp_frac"]]

    @nest.add_target()
    def linearham_final_output(outdir, c):
        linearham_final_output = env.Command(
            [os.path.join(outdir, filename) for filename in
                ["linearham_run.trees", "linearham_run.log", "linearham_run.ess"]],
            [c["linearham_intermediate_output"], c["cluster_fasta_file"]],
            "Rscript --slave --vanilla scripts/run_bootstrap_asr_ess.R" \
                + " $SOURCES" \
                + " " + str(c["linearham_setting"]["burnin_frac"]) \
                + " " + str(c["linearham_setting"]["subsamp_frac"]) \
                + " " + str(1) \
                + " " + str(c["revbayes_setting"]["rng_seed"]) \
                + " $TARGETS")
        env.Depends(linearham_final_output, "scripts/run_bootstrap_asr_ess.R")
        return linearham_final_output

    @nest.add_target()
    def linearham_annotations(outdir, c):
        outbase = os.path.join(outdir, "linearham_annotations")
        linearham_annotations = env.Command(
                [outbase + sfx for sfx in ("_best.yaml", "_all.yaml")],
                [c["cluster_yaml_file"], c["linearham_final_output"][1]],
                "scripts/write_lh_annotations.py $SOURCES --output-base " + outbase)
        env.Depends(linearham_annotations, "scripts/write_lh_annotations.py")
        return linearham_annotations

    @nest.add_target()
    def naive_tabulation(outdir, c):
        outbase = os.path.join(outdir, "aa_naive_seqs")
        naive_tabulation = env.Command(
            outbase + ".fasta", c["linearham_final_output"][0],
            "scripts/tabulate_naive_probs.py $SOURCE --output-base " + outbase)
        env.Depends(naive_tabulation, "scripts/tabulate_naive_probs.py")
        return naive_tabulation

    if options["lineage_unique_ids"] is not None:

        @nest.add_target()
        def lineage_tabulation(outdir, c):
            outbase = os.path.join(outdir, "aa_lineage_seqs")
            lineage_tabulation = env.Command(
                outbase + ".fasta", [c["linearham_final_output"][0], c["naive_tabulation"]],
                "scripts/tabulate_lineage_probs.py $SOURCES" \
                    + " --seed-seq " + c["lineage"]["lineage_seq_unique_id"] \
                    + " --pfilters " + " ".join(str(pfilter) for pfilter in options["asr_pfilters"]) \
                    + " --output-base " + outbase)
            env.Depends(lineage_tabulation, "scripts/tabulate_lineage_probs.py")
            return lineage_tabulation
