"""Some routines for generation workflow"""
from routine import (_remove_mols, _db_check)
from CGRtools.reactor import Reactor
from CGRtools.containers import ReactionContainer
from CGRtools.exceptions import (InvalidAromaticRing,
                                 MappingError)

INITIAL = {"type": "Initial"}
RECONSTRUCTED = {"type": "Reconstructed"}
DECOY = {"type": "Decoy"}


def remove_reagents(reaction):
    """
    Removing unchanging molecules, and checking reaction properties
    :param reaction: input reaction
    :return: ReactionContainer or None in some cases
    """
    try:
        cgr = reaction.compose()
    except ValueError:
        return None
    _remove_mols(reaction, cgr)  # remove unchanging molecules
    try:
        reaction.canonicalize()
        reaction.flush_cache()
    except InvalidAromaticRing:
        return None
    reactants = [x for x in reaction.reactants if x]
    products = [x for x in reaction.products if x]
    if len(reactants) > 0 and len(products) > 0 and _db_check(cgr):  # the num of dyn bonds in a CGR must be > 1
        reaction = ReactionContainer(reactants=reactants,
                                     products=products,
                                     meta=reaction.meta)
        return reaction
    else:
        return None


def not_radical(cgr):
    """
    Checking for charged atoms in a Condensed Graph of Reaction
    :param cgr: Condensed Graph of the input reaction
    :return: bool
    """
    if cgr and cgr.center_atoms:
        if any(x.is_radical or x.p_is_radical for _, x in cgr.atoms()):
            return False
    return True


def containers_split(reaction):
    """
    Separation of MoleculeContainers
    :param reaction: input reaction
    :return: ReactionContainer
    """
    new_reactants = [x for container in reaction.reactants for x in container.split()]
    new_products = [x for container in reaction.products for x in container.split()]
    return ReactionContainer(reactants=new_reactants,
                             products=new_products,
                             meta=reaction.meta)


def generate_reactions(reaction, reactants, rules, max_decoys, limit, doc):
    """
    Accumulation of generated reactions
    :param reaction: input reaction
    :param reactants: reactants of input reaction
    :param rules: rules
    :param max_decoys: max number of reaction to generate
    :param limit: max number of reaction from one transformation
    :param doc: Dict[CGR(reaction), Dict[ReactionContainer, reaction type]}]
    """
    rxn_list = []
    for n, r in enumerate(apply_rules(reactants, rules, limit, max_decoys)):
        if len(rxn_list) == max_decoys:
            break
        new_reaction = ReactionContainer(reactants=reactants,
                                         products=r.products,
                                         meta=r.meta)
        new_reaction.meta.update(reaction.meta)
        try:
            new_reaction.canonicalize()
            new_reaction = remove_reagents(new_reaction)
            if new_reaction is None:
                continue
            new_reaction = containers_split(new_reaction)
            new_reaction.flush_cache()
        except InvalidAromaticRing:
            continue
        try:
            if (new_reaction.compose()).center_atoms:
                try:
                    if str(new_reaction.compose()) not in doc:
                        new_reaction.meta.update(DECOY)
                    elif str(new_reaction.compose()) in doc and \
                            doc[str(new_reaction.compose())]["type"] == "Initial":
                        new_reaction.meta.update(RECONSTRUCTED)
                    else:
                        continue
                except KeyError:
                    continue

                doc.update({str(new_reaction.compose()): {"structure": new_reaction,
                                                          "type": new_reaction.meta["type"]}})
                rxn_list.append(new_reaction)
            else:
                continue
        except MappingError:
            continue


def apply_rules(reactants, rules, limit, max_decoys):
    """
    New reaction generator
    :param reactants: reactants of input reaction
    :param rules: list of rules for input reaction (list[ReactionContainer, ...])
    :param limit: max number of reaction from one transformation
    :param max_decoys: max number of reaction to generate
    :return: yield(ReactionContainer) of generated reaction
    """
    reactors = [Reactor(rule,
                        delete_atoms=True,
                        one_shot=True,
                        automorphism_filter=False) for rule in rules]  # NB! CGRtools v. > 4.1.22
    rxn_list = []
    queue = [r(reactants) for r in reactors]

    while queue:
        if len(rxn_list) >= max_decoys + 10:
            break
        reactor_call = queue.pop(0)
        try:
            rxn_from_apply = []
            for new_reaction in reactor_call:
                if new_reaction in rxn_from_apply or \
                        new_reaction in rxn_list:
                    continue

                rxn_from_apply.append(new_reaction)
                if len(rxn_from_apply) == limit:
                    rxn_list.extend(rxn_from_apply)
                    break
            else:
                rxn_list.extend(rxn_from_apply)
        except KeyError:
            continue
    return rxn_list


def get_rules(reaction):
    """
    Obtaining the rules of reaction transformations
    NB! CGRtools.enumerate_centers() is used
    :param reaction: input reaction
    :return: list[ReactionContainer, ...]
    """
    rules = []

    for n, partial_reaction in enumerate(reaction.enumerate_centers()):  # getting single stages
        if n == 5:
            break  # without an endless loop may appear, cause of many reaction centers

        reactants = []
        products = []

        cleavage = set(partial_reaction.reactants).difference(partial_reaction.products)
        coming = set(partial_reaction.products).difference(partial_reaction.reactants)

        try:
            reaction_center = set(partial_reaction.extended_centers_list[0])  # extended reaction centers
        except IndexError:
            continue

        bare_reaction_center = set(partial_reaction.compose().center_atoms)

        for mol in partial_reaction.reactants:
            group_atoms = reaction_center.intersection(mol)
            if group_atoms:
                group = mol.substructure(group_atoms, as_query=True)
                for i in group_atoms.difference(bare_reaction_center):
                    group._neighbors[i] = ()  # getting rid of the neighbors in rule (NECESSARY)
                reactants.append(group)
        for mol in partial_reaction.products:
            group_atoms = reaction_center.intersection(mol)
            if group_atoms:
                group = mol.substructure(group_atoms, as_query=True)
                for i in group_atoms.difference(bare_reaction_center):
                    group._neighbors[i] = ()  # getting rid of the neighbors in rule (NEC.)
                products.append(group)

        if len(reactants) != 2:
            continue
        rule = ReactionContainer(reactants=reactants,
                                 products=products,
                                 meta=reaction.meta)

        for molecule in rule.molecules():
            molecule._rings_sizes = {x: () for x in molecule._rings_sizes}  # getting rid of the ring sizes info (NEC.)
            molecule._hydrogens = {x: () for x in molecule._hydrogens}  # getting rid of the hydrogen info (NEC.)

        for del_hyb_atom in cleavage:
            for molecule in rule.reactants:
                if del_hyb_atom in molecule:
                    molecule._hybridization[del_hyb_atom] = ()  # getting rid of the hybridization info in react. (NEC.)

        for del_hyb_atom in coming:
            for molecule in rule.products:
                if del_hyb_atom in molecule:
                    molecule._hybridization[del_hyb_atom] = ()  # getting rid of the hybridization info in prod. (NEC.)

        rule.flush_cache()
        rules.append(rule)
    return rules
