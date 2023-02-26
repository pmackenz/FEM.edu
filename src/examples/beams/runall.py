"""
Control file to run all beam examples contained in this source distribution
"""

import __init__

for exmpl in __init__.__all__:
    print(f"   >> {exmpl}")
