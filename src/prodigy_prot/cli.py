"""
Binding affinity predictor based on Intermolecular Contacts (ICs).
"""

import argparse
import logging
import sys
from argparse import RawTextHelpFormatter
from concurrent.futures import ProcessPoolExecutor, as_completed
from io import StringIO
from pathlib import Path

from Bio.PDB.Model import Model

from prodigy_prot.modules.parsers import parse_structure
from prodigy_prot.modules.prodigy import Prodigy

# setup logging
logging.basicConfig(level=logging.INFO, stream=sys.stdout, format="%(message)s")
log = logging.getLogger("Prodigy")


ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ap.add_argument(
    "input_path",
    help="Path to either: \n- Structure in PDB or mmCIF format\n- Directory containing structure files",
)
ap.add_argument(
    "--distance-cutoff",
    type=float,
    default=5.5,
    help="Distance cutoff to calculate ICs",
)
ap.add_argument(
    "--acc-threshold",
    type=float,
    default=0.05,
    help="Accessibility threshold for BSA analysis",
)
ap.add_argument(
    "--temperature",
    type=float,
    default=25.0,
    help="Temperature (C) for Kd prediction",
)
ap.add_argument("--contact_list", action="store_true", help="Output a list of contacts")
ap.add_argument(
    "--pymol_selection",
    action="store_true",
    help="Output a script to highlight the interface (pymol)",
)
ap.add_argument(
    "-q",
    "--quiet",
    action="store_true",
    help="Outputs only the predicted affinity value",
)
ap.add_argument(
    "-np",
    "--number-of-processors",
    type=int,
    action="store",
    help="Number of processors to use (default: 1)",
    default=1,
)
_co_help = """
By default, all intermolecular contacts are taken into consideration,
a molecule being defined as an isolated group of amino acids sharing
a common chain identifier. In specific cases, for example
antibody-antigen complexes, some chains should be considered as a
single molecule.

Use the --selection option to provide collections of chains that should
be considered for the calculation. Separate by a space the chains that
are to be considered _different_ molecules. Use commas to include multiple
chains as part of a single group:

--selection A B => Contacts calculated (only) between chains A and B.
--selection A,B C => Contacts calculated (only) between \
    chains A and C; and B and C.
--selection A B C => Contacts calculated (only) between \
    chains A and B; B and C; and A and C.
"""
sel_opt = ap.add_argument_group("Selection Options", description=_co_help)
sel_opt.add_argument("--selection", nargs="+", metavar=("A B", "A,B C"))


def main():
    args = ap.parse_args()
    log.setLevel(logging.ERROR if args.quiet else logging.INFO)

    struct_path = Path(args.input_path)

    input_list = []
    if struct_path.is_file():
        input_list.append(struct_path)

    elif struct_path.is_dir():
        for input_f in struct_path.glob("*"):
            if Path(input_f).suffix in [".pdb", ".cif", ".ent"]:
                input_list.append(input_f)

    elif not struct_path.exists():
        log.error(f"File {struct_path} does not exist")
        sys.exit(1)

    else:
        log.error(f"Input path {struct_path} is neither a valid file nor a directory")
        sys.exit(1)

    # Collect all tasks
    tasks = []
    for input_f in input_list:
        models, _, _ = parse_structure(str(input_f))
        struct_path = Path(input_f)

        for model in models:
            identifier = f"{struct_path.stem}_model{model.id}"
            tasks.append((model, identifier, args, struct_path))

    # Execute in parallel
    total_tasks = len(tasks)
    if total_tasks == 0:
        log.error("No valid structures found")
        sys.exit(1)
    max_workers = min(args.number_of_processors, total_tasks)
    log.info(f"[+] Executing {total_tasks} task(s) in total")
    if max_workers != args.number_of_processors:
        log.info("[+] Adjusting number of processors based on number of tasks")
        log.info(
            f"[+] Using {max_workers} processor(s) instead of {args.number_of_processors}"
        )

    # Execute and collect results
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_model, *task) for task in tasks]
        for future in as_completed(futures):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                log.error(f"Error processing model: {e}")

    # Sort by identifier, then model.id
    results.sort(key=lambda x: (x[0], x[1]))
    # Print all outputs sequentially
    for identifier, _, output in results:
        print(output, end="")


def process_model(model: Model, identifier: str, args: argparse.Namespace, struct_path):
    """Process a single model"""
    # Capture stdout
    output_buffer = StringIO()

    old_stdout = sys.stdout
    sys.stdout = output_buffer
    try:
        if not args.quiet:
            print("#" * 42)
            print(f"[+] Processing structure {identifier}")
        prodigy = Prodigy(
            model=model,
            name=identifier,
            selection=args.selection,
            temp=args.temperature,
        )
        prodigy.predict(
            distance_cutoff=args.distance_cutoff, acc_threshold=args.acc_threshold
        )
        prodigy.print_prediction(quiet=args.quiet)
    finally:
        sys.stdout = old_stdout

    if args.contact_list:
        contact_list_f = struct_path.with_suffix(".ic")
        prodigy.print_contacts(outfile=str(contact_list_f))

    if args.pymol_selection:
        pymol_script_f = struct_path.with_suffix(".pml")
        prodigy.print_pymol_script(outfile=str(pymol_script_f))

    return identifier, model.id, output_buffer.getvalue()


if __name__ == "__main__":
    sys.exit(main())
