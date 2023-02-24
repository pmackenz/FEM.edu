"""
Main file to run all examples contained in this source distribution
"""
import os, sys
import subprocess

all_examples = [
    'trusses',
    'beams',
    'frames',
    'plates',
    'solids',
    'mixed_structures',
]

if __name__ == "__main__":

    for example_folder in all_examples:
        try:
            os.chdir(example_folder)
        except:
            continue

        print(f"Running examples in: {example_folder}")

        try:
            subprocess.call([sys.executable, "./runall.py"])
        except:
            print("... no runall.py script found in this folder")

        os.chdir('..')
