Search.setIndex({"docnames": ["design", "design/element", "design/material", "design/node", "design/plotter", "design/solver", "design/system", "implementation", "implementation/Element_class", "implementation/FiberMaterial", "implementation/LinearSolver_class", "implementation/LinearTriangle_class", "implementation/Material_class", "implementation/Node_class", "implementation/PlaneStrain", "implementation/Plotter_class", "implementation/Solver_class", "implementation/System_class", "implementation/Truss_class", "index"], "filenames": ["design.rst", "design/element.rst", "design/material.rst", "design/node.rst", "design/plotter.rst", "design/solver.rst", "design/system.rst", "implementation.rst", "implementation/Element_class.rst", "implementation/FiberMaterial.rst", "implementation/LinearSolver_class.rst", "implementation/LinearTriangle_class.rst", "implementation/Material_class.rst", "implementation/Node_class.rst", "implementation/PlaneStrain.rst", "implementation/Plotter_class.rst", "implementation/Solver_class.rst", "implementation/System_class.rst", "implementation/Truss_class.rst", "index.rst"], "titles": ["Program Design", "Element", "Material", "Node", "Plotter", "Solver", "System", "Implementation", "Element classes", "FiberMaterial material class", "Linear Solver class", "LinearTriangle class", "Material classes", "Node class", "PlaneStrain material class", "Plotter class", "Solver class", "System class", "Truss class", "Welcome to the FEM.edu documentation!"], "terms": {"The": [0, 4], "goal": 0, "i": [0, 2, 3, 4, 5, 15, 16], "creat": [0, 4, 5, 6, 15], "an": [0, 1, 5, 6, 16], "object": [0, 1, 4, 5, 6, 15, 17], "orient": 0, "finit": 0, "element": [0, 4, 5, 6, 7, 11, 17, 18, 19], "analysi": 0, "base": [0, 4, 17], "method": 0, "handl": [0, 6], "arbitrari": 0, "truss": [0, 1, 6, 7, 8, 11], "2d": [0, 14], "allow": 0, "load": [0, 3, 5, 6, 17], "can": 0, "plot": [0, 4, 5, 6, 15, 17], "its": [0, 1], "undeform": [0, 4], "deform": [0, 3, 4, 13, 15, 17], "shape": [0, 4], "need": [0, 1, 3, 4, 16], "us": [0, 1, 5, 15], "follow": 0, "node": [0, 1, 4, 5, 6, 7, 17, 19], "materi": [0, 1, 7, 8, 11, 18, 19], "system": [0, 3, 4, 7, 10, 15, 19], "solver": [0, 6, 7, 19], "plotter": [0, 5, 6, 7, 19], "each": [1, 3, 4, 5], "instanc": [1, 3, 5, 6], "repres": [1, 3, 4, 8, 9, 11, 12, 13, 14, 15, 17, 18], "one": [1, 3, 5, 6], "member": 1, "input": [1, 2, 3, 4, 5, 6], "return": [1, 2, 3, 4, 5, 6, 8, 9, 12, 13, 14, 16, 18], "descript": [1, 2, 3, 4, 5, 6], "__init__": [1, 2, 3, 4, 5, 6], "nd0": 1, "nd1": 1, "two": [1, 3, 4], "constructor": [1, 2, 3, 4, 5, 6], "getforc": [1, 5, 8], "list": [1, 3, 4, 5, 6], "1d": [1, 9], "np": [1, 3, 4, 5], "arrai": [1, 3, 4, 5, 15], "A": [1, 2, 9, 10, 12], "nodal": [1, 3, 5, 6, 13, 15], "forc": [1, 3, 4, 5, 6, 15], "shall": [1, 4], "compon": [1, 3], "getstiff": [1, 2, 5, 8, 12], "2": [1, 2, 4, 5], "matrix": 1, "contain": [1, 4, 5], "tangent": [1, 2, 12], "matric": 1, "note": [1, 3], "mai": 1, "have": [1, 3], "chang": 1, "state": [1, 2, 9, 12, 14], "between": 1, "call": 1, "so": [1, 3], "you": 1, "recomput": [1, 5], "everi": 1, "time": 1, "name": [1, 2, 3, 4, 5, 6], "type": [1, 2, 3, 4, 5, 6], "either": 1, "end": [1, 4], "pointer": [1, 6], "comput": [1, 5, 6], "stress": [1, 2, 12, 14], "modulu": [1, 2], "hold": [1, 2, 3, 5, 6], "p0": 1, "p1": 1, "kt": 1, "stiff": [1, 2, 5, 12], "all": [1, 2, 5, 6], "equat": [1, 2, 3, 5], "bf": [1, 3, 5, 10, 13], "l": 1, "x": [1, 3, 4, 13], "_1": [1, 3], "_0": [1, 3], "ell": 1, "n": [1, 9, 14], "frac": [1, 2], "1": [1, 2, 3, 5, 6, 9, 12, 13, 14, 17], "strain": [1, 2, 9, 12, 14], "varepsilon": [1, 2], "cdot": 1, "u": [1, 3, 4, 5, 6, 10, 13], "f": [1, 2, 5, 6, 13], "sigma": [1, 2, 12], "setstrain": [1, 2, 9, 12, 14], "ep": [1, 2, 9, 12, 14], "getstress": [1, 2, 12], "vector": [1, 3, 4, 5, 13, 15, 17], "p": [1, 5, 6, 10], "e": [1, 2, 9, 12, 14], "k": [1, 5, 10], "e_t": [1, 2], "otim": 1, "find": 1, "_": [1, 5], "00": 1, "11": [1, 5], "01": 1, "10": [1, 2], "thi": [2, 3, 4, 5, 10, 16], "provid": [2, 3, 9, 12, 14], "demonstr": 2, "exampl": 2, "paramet": [2, 9, 12, 13, 14, 15, 17], "0": [2, 3, 5, 6, 9, 12, 13, 14, 17], "set": [2, 3, 4, 5], "initi": [2, 3, 4, 13, 16], "intern": [2, 3, 15], "getarea": 2, "cross": 2, "section": 2, "area": 2, "from": [2, 5], "request": [2, 5, 6, 12], "axial": [2, 9, 12, 14], "updat": [2, 5, 6, 9, 12, 14], "user": [2, 4, 9, 12, 14], "valu": [2, 4, 9, 12, 14], "param": [2, 9, 12, 14], "dict": 2, "default": [2, 3, 4], "100": 2, "nu": [2, 9, 12, 14], "fy": [2, 9, 12, 14], "0e30": 2, "moe": 2, "poisson": 2, "": [2, 5], "ratio": 2, "yield": 2, "plastic_strain": 2, "float": [2, 3], "sig": 2, "current": [2, 3], "et": [2, 12], "materil": 2, "elast": 2, "trial": 2, "varepsilon_p": 2, "check": 2, "f_y": 2, "IF": 2, "ge": 2, "3": 2, "delta": 2, "text": [2, 17], "sign": 2, "y": [3, 4, 13], "coordin": [3, 4], "point": [3, 4], "posit": [3, 4, 5, 6, 13], "displac": [3, 4, 5, 6, 13, 17], "zero": [3, 5], "fixdof": [3, 13], "idx": [3, 13], "degre": 3, "freedom": 3, "dof": 3, "flag": 3, "accordingli": 3, "isfix": [3, 13], "true": [3, 15], "fals": [3, 4, 15], "test": 3, "function": [3, 16], "fix": [3, 5], "otherwis": 3, "setdisp": [3, 5, 13], "v": [3, 4, 5, 13], "overwrit": 3, "getdisp": [3, 13], "getpo": [3, 13], "getdeformedpo": [3, 13], "factor": [3, 5, 6, 13, 17], "magnifi": 3, "would": 3, "good": 3, "none": [3, 4, 15], "given": [3, 4, 15], "addload": 3, "px": 3, "py": 3, "add": [3, 4, 5, 6, 15], "setload": 3, "replac": [3, 4], "getload": [3, 5], "po": 3, "index": [3, 5, 19], "int": 3, "addnod": [3, 5, 6, 17], "thisnod": 3, "team": 3, "disp": [3, 4, 5, 15], "fixiti": 3, "th": 3, "sensibl": 4, "setmesh": [4, 15], "vert": [4, 15], "line": [4, 15], "indic": [4, 5, 15], "self": 4, "vertic": 4, "inform": 4, "setdisplac": [4, 15], "setvalu": [4, 15], "val": [4, 15], "displacementplot": [4, 15], "file": [4, 15], "string": 4, "show": 4, "black": 4, "model": [4, 5, 6, 17], "red": 4, "If": [4, 5, 15], "save": 4, "copi": 4, "valueplot": [4, 15], "color": [4, 15], "colormap": 4, "colorbar": 4, "legend": 4, "pair": 4, "start": 4, "respect": 4, "must": 4, "ident": 4, "entri": 4, "newnod": [5, 6, 17], "your": [5, 6], "addel": [5, 6, 17], "newelem": [5, 6], "solv": [5, 6, 16, 17], "assembl": [5, 6, 16], "k_t": [5, 6], "loop": [5, 6], "through": [5, 6], "unbalanc": [5, 6], "r": [5, 6, 15], "collect": [5, 6], "info": [5, 6], "send": [5, 6], "report": [5, 6, 17], "print": [5, 6, 17], "summari": [5, 6, 17], "size": 5, "ha": 5, "node0": [5, 8, 11, 18], "elem": 5, "node1": [5, 8, 11, 18], "j": 5, "local": 5, "d": 5, "o": 5, "belong": 5, "global": 5, "4": 5, "m": 5, "assembli": 5, "sy": 5, "over": 5, "should": 5, "0e20": 5, "step": 5, "5": 5, "6": 5, "do": 5, "repeat": 5, "7": 5, "everyth": 5, "wa": 5, "done": 5, "correctli": 5, "support": 5, "reaction": 5, "free": 5, "numer": 5, "class": [7, 19], "deriv": 7, "lineartriangl": [7, 8], "fibermateri": [7, 12], "planestrain": [7, 12], "linear": [7, 16], "src": [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], "abstract": [8, 12, 16], "singl": [8, 11, 13, 18], "gener": [8, 12, 16], "getaxialforc": [8, 18], "updatest": [8, 11, 18], "1e": [9, 12, 14], "30": [9, 12, 14], "fiber": 9, "tensor": [9, 12, 14], "linearsolv": 10, "impli": 10, "node2": 11, "magnif": [13, 17], "t": 14, "plane": 14, "addforc": 15, "ax": 15, "shown": 15, "axi": 15, "which": 15, "store": 15, "proper": 15, "extens": 15, "desir": 15, "format": 15, "png": 15, "pdf": 15, "filenam": 15, "str": 15, "setreact": 15, "identifi": 15, "magnitud": 15, "defin": 16, "interfac": 16, "ani": 16, "implement": [16, 19], "describ": 16, "reset": [16, 17], "newel": 17, "resetal": 17, "resetdisp": 17, "resetload": 17, "program": 19, "design": 19, "modul": 19, "search": 19, "page": 19}, "objects": {"src": [[13, 0, 0, "-", "Node"], [15, 0, 0, "-", "Plotter"], [17, 0, 0, "-", "System"]], "src.Node": [[13, 1, 1, "", "Node"]], "src.Node.Node": [[13, 2, 1, "", "fixDOF"], [13, 2, 1, "", "getDeformedPos"], [13, 2, 1, "", "getDisp"], [13, 2, 1, "", "getPos"], [13, 2, 1, "", "isFixed"], [13, 2, 1, "", "setDisp"]], "src.Plotter": [[15, 1, 1, "", "Plotter"]], "src.Plotter.Plotter": [[15, 2, 1, "", "addForces"], [15, 2, 1, "", "displacementPlot"], [15, 2, 1, "", "setDisplacements"], [15, 2, 1, "", "setMesh"], [15, 2, 1, "", "setReactions"], [15, 2, 1, "", "setValues"], [15, 2, 1, "", "valuePlot"]], "src.System": [[17, 1, 1, "", "System"]], "src.System.System": [[17, 2, 1, "", "addElement"], [17, 2, 1, "", "addNode"], [17, 2, 1, "", "plot"], [17, 2, 1, "", "report"], [17, 2, 1, "", "resetAll"], [17, 2, 1, "", "resetDisp"], [17, 2, 1, "", "resetLoad"], [17, 2, 1, "", "solve"]], "src.elements": [[8, 0, 0, "-", "Element"], [11, 0, 0, "-", "LinearTriangle"], [18, 0, 0, "-", "Truss"]], "src.elements.Element": [[8, 1, 1, "", "Element"]], "src.elements.Element.Element": [[8, 2, 1, "", "getAxialForce"], [8, 2, 1, "", "getForce"], [8, 2, 1, "", "getStiffness"], [8, 2, 1, "", "updateState"]], "src.elements.LinearTriangle": [[11, 1, 1, "", "LinearTriangle"]], "src.elements.LinearTriangle.LinearTriangle": [[11, 2, 1, "", "updateState"]], "src.elements.Truss": [[18, 1, 1, "", "Truss"]], "src.elements.Truss.Truss": [[18, 2, 1, "", "getAxialForce"], [18, 2, 1, "", "updateState"]], "src.materials": [[9, 0, 0, "-", "FiberMaterial"], [12, 0, 0, "-", "Material"], [14, 0, 0, "-", "PlaneStrain"]], "src.materials.FiberMaterial": [[9, 1, 1, "", "FiberMaterial"]], "src.materials.FiberMaterial.FiberMaterial": [[9, 2, 1, "", "setStrain"]], "src.materials.Material": [[12, 1, 1, "", "Material"]], "src.materials.Material.Material": [[12, 2, 1, "", "getStiffness"], [12, 2, 1, "", "getStress"], [12, 2, 1, "", "setStrain"]], "src.materials.PlaneStrain": [[14, 1, 1, "", "PlaneStrain"]], "src.materials.PlaneStrain.PlaneStrain": [[14, 2, 1, "", "setStrain"]], "src.solver": [[10, 0, 0, "-", "LinearSolver"], [16, 0, 0, "-", "Solver"]], "src.solver.LinearSolver": [[10, 1, 1, "", "LinearSolver"]], "src.solver.Solver": [[16, 1, 1, "", "Solver"]], "src.solver.Solver.Solver": [[16, 2, 1, "", "assemble"], [16, 2, 1, "", "initialize"], [16, 2, 1, "", "reset"], [16, 2, 1, "", "solve"]]}, "objtypes": {"0": "py:module", "1": "py:class", "2": "py:method"}, "objnames": {"0": ["py", "module", "Python module"], "1": ["py", "class", "Python class"], "2": ["py", "method", "Python method"]}, "titleterms": {"program": 0, "design": 0, "element": [1, 2, 8], "class": [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], "method": [1, 2, 3, 4, 5, 6], "variabl": [1, 2, 3, 4, 5, 6], "materi": [2, 9, 12, 14], "node": [3, 13], "plotter": [4, 15], "solver": [5, 10, 16], "system": [5, 6, 17], "implement": 7, "deriv": [8, 12, 16], "fibermateri": 9, "parent": [9, 10, 11, 14, 18], "doc": [9, 10, 11, 14, 18], "linear": 10, "lineartriangl": 11, "planestrain": 14, "truss": 18, "welcom": 19, "fem": 19, "edu": 19, "document": 19, "content": 19, "indic": 19, "tabl": 19}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx": 56}})