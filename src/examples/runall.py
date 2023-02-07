"""
Main file to run all examples contained in this source distribution
"""
import os, sys

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

        print(f"entering {example_folder}")

        os.chdir('..')
