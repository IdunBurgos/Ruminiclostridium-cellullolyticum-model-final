import logging
import re
from functools import partial
from types import ModuleType
from typing import TYPE_CHECKING, Dict, List, NamedTuple, Optional, Tuple, Union
from warnings import warn

import optlang
import pandas as pd
from optlang.interface import (
    FEASIBLE,
    INFEASIBLE,
    ITERATION_LIMIT,
    NUMERIC,
    OPTIMAL,
    SUBOPTIMAL,
    TIME_LIMIT,
)
from optlang.symbolics import Basic, Zero

from cobra.exceptions import (
    OPTLANG_TO_EXCEPTIONS_DICT,
    OptimizationError,
    SolverNotFound,
)
from cobra.util.context import get_context

from cobra.util import fix_objective_as_constraint

def add_lexicographic_constraints(
    model: "Model",
    objectives: List["Reaction"],
    objective_direction: Union[str, List[str]] = "max",
) -> pd.Series:
    """Successively optimize separate targets in a specific order.

    For each objective, optimize the model and set the optimal value as a
    constraint. Proceed in the order of the objectives given. Due to the
    specific order this is called lexicographic FBA [1]_. This procedure
    is useful for returning unique solutions for a set of important
    fluxes. Typically this is applied to exchange fluxes.

    Parameters
    ----------
    model : cobra.Model
        The model to be optimized.
    objectives : list of cobra.Reaction
        A list of reactions (or objectives) in the model for which unique
        fluxes are to be determined.
    objective_direction : str or list of str, optional
        The desired objective direction for each reaction (if a list) or
        the objective direction to use for all reactions (default "max").

    Returns
    -------
    pandas.Series
        A pandas Series containing the optimized fluxes for each of the
        given reactions in `objectives`.

    References
    ----------
    .. [1] Gomez, Jose A., Kai Höffner, and Paul I. Barton.
    “DFBAlab: A Fast and Reliable MATLAB Code for Dynamic Flux Balance
    Analysis.” BMC Bioinformatics 15, no. 1 (December 18, 2014): 409.
    https://doi.org/10.1186/s12859-014-0409-8.

    """

    if type(objective_direction) is not list:
        objective_direction = [objective_direction] * len(objectives)

    constraints = []
    for rxn_id, obj_dir in zip(objectives, objective_direction):
        print(rxn_id)
        model.objective = model.reactions.get_by_id(rxn_id)
        model.objective_direction = obj_dir
        constraints.append(fix_objective_as_constraint(model))

    return pd.Series(constraints, index=objectives)