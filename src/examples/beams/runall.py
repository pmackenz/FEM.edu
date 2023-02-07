"""
Control file to run all beam examples contained in this source distribution
"""
import os, sys

all_examples = [
    '',
]

if __name__ == "__main__":
    for example_folder in all_examples:
        try:
            os.chdir(example_folder)
        except:
            continue

        print(f"entering {example_folder}")

        os.chdir('..')
