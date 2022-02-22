"""Decoy generation workflow"""

# Import relevant packages
from CGRtools.files import RDFRead, RDFWrite
from CGRtools.exceptions import *
from ..util.utils import (generate_reactions, remove_reagents,
                          containers_split, not_radical,
                          get_rules)
from ..util.routine import (_save_log, _util_file)
from datetime import date

import os
import math
import multiprocessing
import pickle
import time
import argparse


def main(n_proc):
    """
    Main generation routine

    =======================   ===================================
    Reactions in input file   1 316 541
    Data                      reactions, type: ReactionContainer
    =======================   ===================================

    Parameters
    ----------
    :param n_proc: Num pool worker

    Returns
    -------
    :return: None
    """
    path = "{}/data/decoyGeneration/".format(os.path.abspath(os.path.join(os.getcwd(), os.pardir)))
    with open("{}Config.pickle".format(path), "rb") as configFile:
        config_list = pickle.load(configFile)
        name_in, templates_fn, name_out, batch, max_decoys, limit, v, log = config_list

    with open(templates_fn, "rb") as pkl:
        templates = pickle.load(pkl)

    with RDFRead(name_in, indexable=True) as data:

        name = "Worker-{}".format(n_proc)
        id_start = n_proc * batch
        id_stop = id_start + batch
        temp = []

        # LOGGING
        if log:
            log_filename = "{}GENERATE_DECOYS_LOG.txt".format(path)
        if log:
            _save_log(log_filename, str("Process {} initialized\n".format(name)))
        start_time = time.time()

        for n, reaction in enumerate(data[id_start: id_stop], start=1):
            doc = {}
            reaction = remove_reagents(reaction)
            if reaction is not None:
                reaction = containers_split(reaction)
                if len(reaction.reactants) == 2:
                    cgr = reaction.compose()
                    if not_radical(cgr):
                        reaction.meta.update(INITIAL)
                        doc.update({str(reaction.compose()): {"structure": reaction,
                                                              "type": reaction.meta["type"]}})
                        rules = get_rules(reaction)
                        if rules:
                            generate_reactions(reaction, reaction.reactants, rules,
                                               max_decoys, limit, doc)
                        else:
                            if log:
                                _save_log(log_filename,
                                          str("Failed to get strict templates for reaction with ID {}\n".format(
                                              reaction.meta["Reaction_ID"])))

                        generate_reactions(reaction, reaction.reactants, templates,
                                           max_decoys, limit, doc)
                        if any([True if v["structure"].meta["type"].startswith("Initial") else False for v in
                                doc.values()]):
                            if v:
                                print("Reaction with ID {} was not recovered\n".format(reaction.meta["Reaction_ID"]))
                            if log:
                                _save_log(log_filename,
                                          str("Reaction with ID: {} was not recovered\n".format(
                                              reaction.meta["Reaction_ID"])))
                            continue

                        if any([True if v["structure"].meta["type"].startswith("Reconstructed") else False for v in
                                doc.values()]):
                            if v:
                                print("Reaction with ID {} was successfully recovered\n".format(
                                        reaction.meta["Reaction_ID"]))
                            if log:
                                _save_log(log_filename,
                                          str("Reaction with ID: {} was successfully recovered\n".format(
                                              reaction.meta["Reaction_ID"])))
                            temp.extend([x["structure"] for x in doc.values()])
            if batch < 10000:
                continue
            if n % (batch / 10) == 0 and n != batch:
                iter_time = time.time()
                if v:
                    print("Process {} stepped over: {} by {}s\n".format(name, n, iter_time - start_time))
                if log:
                    _save_log(log_filename,
                              str("Process {} stepped over: {} by {}s\n".format(name, n, iter_time - start_time)))
                with open("{}/data/decoyGeneration/{}".format(os.path.join(os.getcwd(), os.pardir),
                                                              name_out), "a") as w, RDFWrite(w) as rdf:
                    for item in temp:
                        rdf.write(item)
                temp = []
        else:
            end_time = time.time()
            if v:
                print("Process {} finished batch processing in time: {}s\n".format(name, end_time - start_time))
            if log:
                _save_log(log_filename,
                          str("Process {} finished batch processing in time: {}s\n".format(name,
                                                                                           end_time - start_time)))
            with open("{}/data/decoyGeneration/{}".format(os.path.join(os.getcwd(), os.pardir),
                                                          name_out), "a") as w, RDFWrite(w) as rdf:
                for item in temp:
                    rdf.write(item)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('name_in', type=str,
                        help='Path to the main rdf file')
    parser.add_argument('template_pkl', type=str,
                        help='Path to the templates pickle file')
    parser.add_argument('-num_proc', type=int, default=1,
                        help='Number of processor cores to be used; defaults to 1')
    parser.add_argument('-name_out', type=str, default='Decoys from {}.rdf'.format(date.today()),
                        help='Output rdf file name; defaults to: Decoys from {current date}')
    parser.add_argument('-batch', type=int, default=10000,
                        help='Batch size(the number of studied reactions per pass before recording); defaults to 10000')
    parser.add_argument('-v', type=bool, default=False,
                        help='Verbose printing; defaults to False')
    parser.add_argument('-n', '--num', type=int, default=50,
                        help='The maximum number of decoys that can be created; defaults to 50')
    parser.add_argument('-l', '--lim', type=int, default=5,
                        help='The maximum number of template applying; defaults to 5')
    parser.add_argument('--count', type=int, default=1000,
                        help='Number of templates to be used; defaults to 1000')
    parser.add_argument('--log', type=bool, default=True,
                        help='Whether to log wall times / errors check / etc., default True')
    args = parser.parse_args()
    log = bool(args.log)

    path = "{}/data/decoyGeneration/".format(os.path.abspath(os.path.join(os.getcwd(), os.pardir)))

    INITIAL = {"type": "Initial"}
    RECONSTRUCTED = {"type": "Reconstructed"}
    DECOY = {"type": "Decoy"}

    _util_file(str("{}{}".format(path, args.name_out)))  # delete the rdf file if it was created earlier
    if log:
        _util_file("{}GENERATE_DECOYS_LOG.txt".format(path))

    with open("{}Config.pickle".format(path), "wb") as config:
        config_list = [
            str(args.name_in),
            str(args.template_pkl),
            str(args.name_out),
            int(args.batch),
            int(args.num),
            int(args.lim),
            bool(args.v),
            bool(args.log)
        ]
        pickle.dump(config_list, config)

    with RDFRead(str(args.name_in), indexable=True) as file:
        chunk = math.ceil(len(file) / 1000)

    with multiprocessing.Pool(processes=int(args.num_proc)) as pool:
        pool.map(main, list(range(chunk)))
