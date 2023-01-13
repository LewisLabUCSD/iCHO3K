from warnings import warn
from six import string_types
import re

import cobra

# This regular expression finds any single letter compartment enclosed in
# square brackets at the beginning of the string. For example [c] : foo --> bar
compartment_finder = re.compile("^\s*(\[[A-Za-z]\])\s*:*")
def set_stoichiometry_from_string(reaction, reaction_str, verbose=True):
    if reaction._model is None:
        warn("no model found")
        model = None
    else:
        model = reaction._model
    original_str = "" + reaction_str # copy
    found_compartments = compartment_finder.findall(reaction_str)
    if len(found_compartments) == 1:
        compartment = found_compartments[0]
        reaction_str = compartment_finder.sub("", reaction_str)
    else:
        compartment = ""

    
    if "<->" in reaction_str:
        reaction.lower_bound = -1000
        reactant_str, product_str = reaction_str.split("<->")
    elif "<==>" in reaction_str:
        reaction.lower_bound = -1000
        reactant_str, product_str = reaction_str.split("<==>")
    elif "-->" in reaction_str:
        reactant_str, product_str = reaction_str.split("-->")
    elif "->" in reaction_str:
        reactant_str, product_str = reaction_str.split("->")
    elif "<--" in reaction_str:
        reaction.upper_bound = 0
        reaction.lower_bound = -1000
        reactant_str, product_str = reaction_str.split("<--")
    elif "<-" in reaction_str:
        reaction.upper_bound = 0
        reaction.lower_bound = -1000
        reactant_str, product_str = reaction_str.split("<-")
    else:
        raise ValueError("no suitable arrow found in '%s'" % reaction_str)

    for substr, factor in ((reactant_str, -1), (product_str, 1)):
        substr = substr.strip()
        if len(substr) == 0:
            continue
        for term in substr.split("+"):
            term = term.strip()
            if term.lower() == "nothing":
                continue
            if " " in term:
                num_str, met_id = term.split()
                num = float(num_str.lstrip("(").rstrip(")")) * factor
            else:
                met_id = term
                num = factor
            met_id += compartment
            try:
                met = model.metabolites.get_by_id(met_id)
            except KeyError:
                if verbose:
                    print("unknown metabolite '%s' created" % met_id)
                met = cobra.Metabolite(met_id)
            reaction.add_metabolites({met: num})