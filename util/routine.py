"""Routine functions"""
from CGRtools.files import (RDFRead, RDFWrite)
from tqdm import tqdm

import pickle
import os


def _load_pkl(filename, n):
    with open("{}_{}.pickle".format(filename, n), "rb") as f:
        for i in range(2 ** 64):
            try:
                yield pickle.load(f)
            except EOFError:
                break


def _dump_pkl(filename, to_save, n):
    with open("{}_{}.pickle".format(filename, n), "ab") as f:
        pickle.dump(to_save, f)
    return


def _db_check(cgr):
    if len(cgr.center_bonds) > 1 and len(cgr.center_bonds) != 0:
        return True
    else:
        return False


def _remove_mols(reaction, cgr):
    to_remove = []
    for n, i in enumerate(cgr.split()):
        if any([z[1].charge != z[1].p_charge for z in i.atoms()]) or \
                any([z[2].order != z[2].p_order for z in i.bonds()]):
            continue
        to_remove.extend([x for x, _ in i.atoms()])
    for mol in reaction.molecules():
        for z in to_remove:
            try:
                mol.delete_atom(z)
            except KeyError:
                continue


def _util_file(filename):
    try:
        os.remove(filename)
    except FileNotFoundError:
        pass


def _save_log(filename, string):
    with open(filename, "a") as flog:
        flog.write(string)


def RDFclean(RDFfilename, log, dump_size, dump_fn, v):
    """
    Dump pickling and removing duplicates
    :param RDFfilename: RDF file to check duplicates and examine
    :param log: log if necessary
    :param dump_size: size of dict to dump
    :param dump_fn: outer file name
    :param v: printing if necessary
    """
    with RDFRead(RDFfilename, indexable=True) as file:

        to_save = dict()
        log_filename = "CLEANING_LOG.txt"
        flag = True
        flag_pass = True

        for n, reaction in enumerate(tqdm(file), start=1):
            if n != 2210001 and flag_pass:
                continue
            else:
                flag_pass = False
            if flag:
                num = n - 1
                flag = False
            try:
                if str(reaction.compose()) not in to_save:
                    to_save.update({str(reaction.compose()): reaction})
                else:
                    if reaction.meta["type"].startswith("Reconstructed"):
                        if v: print("Replaced by reconstructed: {}".format(reaction.meta["Reaction_ID"]))
                        if log: _save_log(log_filename,
                                          str("Replaced by reconstructed: {}\n".format(reaction.meta["Reaction_ID"])))
                        del to_save[str(reaction.compose())]
                        to_save.update({str(reaction.compose()): reaction})
                    elif reaction.meta["type"].startswith("Decoy"):
                        if v: print("Founded duplicate decoy: {}, real id: {}".format(reaction.meta["Reaction_ID"],
                                                                                      to_save[
                                                                                          str(reaction.compose())].meta[
                                                                                          "Reaction_ID"]))
                        if log: _save_log(log_filename, str("Founded duplicate decoy: {}, real id: {}\n".format(
                            reaction.meta["Reaction_ID"],
                            to_save[str(reaction.compose())].meta["Reaction_ID"])))
                        continue
            except Exception as e:
                if v: print("{} was occurred, number: {}, rxn_ID: {}\n".format(e, n, reaction.meta["Reaction_ID"]))
                if log: _save_log(log_filename, str("{} was occurred, number: {}, rxn_ID: {}\n".format(e, n,
                                                                                                      reaction.meta[
                                                                                                          "Reaction_ID"]
                                                                                                       )))
                continue
            if n % dump_size == 0:
                flag = True
                try:
                    for rxn_dict in _load_pkl(dump_fn, num):
                        for duplicate in set(to_save).intersection(set(rxn_dict)):
                            try:
                                if to_save[duplicate].meta["type"].startswith("Reconstructed"):
                                    del rxn_dict[duplicate]
                                elif to_save[duplicate].meta["type"].startswith("Decoy"):
                                    del to_save[duplicate]
                            except KeyError:
                                del to_save[duplicate]
                                continue
                        _dump_pkl(dump_fn, rxn_dict, n)
                    else:
                        _dump_pkl(dump_fn, to_save, n)
                        to_save = dict()
                        _util_file("{}_{}.pickle".format(dump_fn, num))
                except (FileNotFoundError, UnboundLocalError):
                    _dump_pkl(dump_fn, to_save, n)
                    to_save = dict()
        else:
            for rxn_dict in _load_pkl(dump_fn, num):
                for duplicate in set(to_save).intersection(set(rxn_dict)):
                    try:
                        if to_save[duplicate].meta["type"].startswith("Reconstructed"):
                            del rxn_dict[duplicate]
                        elif to_save[duplicate].meta["type"].startswith("Decoy"):
                            del to_save[duplicate]
                    except KeyError:
                        del to_save[duplicate]
                        continue
                _dump_pkl(dump_fn, rxn_dict, n)
            else:
                _dump_pkl(dump_fn, to_save, n)
                _util_file("{}_{}.pickle".format(dump_fn, num))

def Compile(input_file, output_file):
    with open(input_file, "rb") as file, \
            open(output_file, "ab") as new_file, \
            open("Nonrecon.pickle", "ab") as Nonrec_file:
        counter = {}
        len_list = []
        nonrec_counter = {}
        num_examined = 0
        id_len = 0
        init = "Not found yet"

        for i in range(2 ** 64):
            try:
                data = pickle.load(file)
                len_list.append(len(data))
                if i == 0:
                    pass
                else:
                    id_len += 1
            except EOFError:
                break

            for n, reaction in enumerate(tqdm(data.values())):
                num_examined += 1
                try:
                    if reaction.meta["Reaction_ID"] != init:
                        try:
                            if counter[init]["Recon_str"] != 0:
                                pass
                            else:
                                nonrec_counter.update({init: counter[init]})
                                del counter[init]
                        except KeyError:
                            print("KeyError occurred: {}".format(reaction.meta["Reaction_ID"]))
                            pass

                        if num_examined >= len_list[id_len - 1]:
                            pickle.dump(counter, new_file)
                            pickle.dump(nonrec_counter, Nonrec_file)
                            num_examined = 0
                            counter = {}
                            nonrec_counter = {}

                        init = reaction.meta["Reaction_ID"]
                        num = 0
                        num_random = 0
                        num_strict = 0
                        structures = []
                        recon_str = 0
                        rand_ids = []
                    structures.append(reaction)
                    num += 1
                    if not reaction.meta["type"].startswith("Reconstructed"):
                        if "Rule_ID" in reaction.meta:
                            num_random += 1
                            rand_ids.append(reaction.meta["Rule_ID"])
                        else:
                            num_strict += 1
                    else:
                        recon_str = reaction
                    counter.update({reaction.meta["Reaction_ID"]:
                        {
                            "Recon_str": recon_str,
                            "Numbers": num,
                            "Structures": structures,
                            "Number of decoys from strict": num_strict,
                            "Number of decoys from random": num_random,
                            "Random rules ID's": rand_ids
                        }
                    })
                except Exception as e:
                    print("{} was occurred, number: {}".format(e, n))
                    continue
