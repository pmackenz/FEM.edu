.
├── Contents.md
├── README.md
├── sphinx
│  ├── build
│  ├── docs
│  ├── make.bat
│  ├── Makefile
│  └── source
│     ├── conf.py
│     ├── design
│     │  ├── element.rst
│     │  ├── material.rst
│     │  ├── node.rst
│     │  ├── plotter.rst
│     │  ├── solver.rst
│     │  └── system.rst
│     ├── design.rst
│     ├── examples
│     │  ├── beam_examples.rst
│     │  ├── beams
│     │  │  ├── beam01
│     │  │  │  ├── beam01.rst
│     │  │  │  ├── beam01_deformed.png
│     │  │  │  ├── beam01_moment.png
│     │  │  │  └── beam01_shear.png
│     │  │  └── beam02
│     │  │     ├── beam02.rst
│     │  │     ├── beam02_deformed.png
│     │  │     ├── beam02_moment.png
│     │  │     └── beam02_shear.png
│     │  ├── continuum_examples.rst
│     │  ├── frame_examples.rst
│     │  ├── frames
│     │  │  ├── frame01
│     │  │  │  ├── frame01.rst
│     │  │  │  ├── frame1_buckling_mode0.png
│     │  │  │  ├── frame1_deformed.png
│     │  │  │  ├── frame1_force.png
│     │  │  │  ├── frame1_moment.png
│     │  │  │  └── frame1_shear.png
│     │  │  ├── frame02
│     │  │  │  ├── frame02.rst
│     │  │  │  ├── frame2_buckling_mode0.png
│     │  │  │  ├── frame2_deformed.png
│     │  │  │  ├── frame2_force.png
│     │  │  │  ├── frame2_moment.png
│     │  │  │  └── frame2_shear.png
│     │  │  ├── frame03
│     │  │  │  ├── frame03.rst
│     │  │  │  ├── frame3_buckling_mode0.png
│     │  │  │  ├── frame3_deformed.png
│     │  │  │  ├── frame3_force.png
│     │  │  │  ├── frame3_moment.png
│     │  │  │  ├── frame3_shear.png
│     │  │  │  └── frame3_stability_check.png
│     │  │  └── frame04
│     │  │     ├── frame04.rst
│     │  │     ├── frame04_buckling_mode0.png
│     │  │     ├── frame04_buckling_mode1.png
│     │  │     ├── frame04_buckling_mode2.png
│     │  │     ├── frame04_buckling_mode3.png
│     │  │     ├── frame04_deformed.png
│     │  │     ├── frame04_force.png
│     │  │     ├── frame04_moment.png
│     │  │     ├── frame04_shear.png
│     │  │     └── frame4_stability_check.png
│     │  ├── mixed_structures_examples.rst
│     │  ├── plate_examples.rst
│     │  ├── truss_examples.rst
│     │  └── trusses
│     │     ├── truss01
│     │     │  ├── truss01.rst
│     │     │  ├── truss01_deformed.png
│     │     │  ├── truss01_deformed_b.png
│     │     │  ├── truss01_forces.png
│     │     │  └── truss01_forces_b.png
│     │     ├── truss02
│     │     │  ├── truss02.rst
│     │     │  ├── truss02_deformed.png
│     │     │  └── truss02_forces.png
│     │     ├── truss03
│     │     │  ├── truss03.rst
│     │     │  ├── truss03_deformed_a.png
│     │     │  └── truss03_deformed_b.png
│     │     └── truss04
│     │        └── truss04.rst
│     ├── examples.rst
│     ├── general
│     │  ├── about.rst
│     │  └── welcome.rst
│     ├── images
│     │  ├── FEMedu-CoverPicture.png
│     │  ├── FEMedu-CoverPicture.psd
│     │  ├── galik_william.jpg
│     │  ├── jordan-seawright.jpg
│     │  ├── peter-mackenzie-helnwein.jpeg
│     │  └── tatsu-sweet.jpeg
│     ├── implementation
│     │  ├── Beam2D_class.rst
│     │  ├── DOF_class.rst
│     │  ├── ElasticSection.rst
│     │  ├── Element_class.rst
│     │  ├── ElementPlotter3D_class.rst
│     │  ├── ElementPlotter_class.rst
│     │  ├── FiberMaterial.rst
│     │  ├── FiberSection.rst
│     │  ├── Frame2D_class.rst
│     │  ├── LinearSolver_class.rst
│     │  ├── LinearTriangle_class.rst
│     │  ├── Material_class.rst
│     │  ├── NewtonRaphsonSolver_class.rst
│     │  ├── Node_class.rst
│     │  ├── PlaneStrain.rst
│     │  ├── PlaneStress.rst
│     │  ├── PlateSection.rst
│     │  ├── Plot_Support_Classes.rst
│     │  ├── Plotter3D_class.rst
│     │  ├── Plotter_class.rst
│     │  ├── Sections.rst
│     │  ├── Solver_class.rst
│     │  ├── System_class.rst
│     │  ├── Transformation_class.rst
│     │  └── Truss_class.rst
│     ├── implementation.rst
│     └── index.rst
└── src
   ├── demo.py
   ├── Demo_notebooks
   │  ├── FEM beam example 01.ipynb
   │  ├── FEM beam example 02.ipynb
   │  ├── FEM frame example 01.ipynb
   │  ├── FEM frame example 01b.ipynb
   │  ├── FEM frame example 01c.ipynb
   │  ├── FEM frame example 02.ipynb
   │  ├── FEM frame example 03.ipynb
   │  ├── FEM frame example 04.ipynb
   │  ├── FEM frame example 05.ipynb
   │  ├── FEM truss example 01.ipynb
   │  ├── FEM truss example 02.ipynb
   │  ├── FEM truss example 03.ipynb
   │  └── Theory.ipynb
   ├── femedu
   │  ├── __init__.py
   │  ├── domain
   │  │  ├── __init__.py
   │  │  ├── DoF.py
   │  │  ├── Node.py
   │  │  ├── System.py
   │  │  └── Transformation.py
   │  ├── elements
   │  │  ├── __init__.py
   │  │  ├── Beam2D.py
   │  │  ├── DrawElement.py
   │  │  ├── Element.py
   │  │  ├── Frame2D.py
   │  │  ├── LinearTriangle.py
   │  │  └── Truss.py
   │  ├── examples
   │  │  ├── __init__.py
   │  │  ├── beams
   │  │  │  ├── __init__.py
   │  │  │  ├── beam01.py
   │  │  │  ├── beam02.py
   │  │  │  └── runall.py
   │  │  ├── Example.py
   │  │  ├── frames
   │  │  │  ├── __init__.py
   │  │  │  ├── frame01.py
   │  │  │  ├── frame02.py
   │  │  │  ├── frame03.py
   │  │  │  ├── frame04.py
   │  │  │  └── runall.py
   │  │  ├── mixed
   │  │  │  ├── __init__.py
   │  │  │  └── runall.py
   │  │  ├── plates
   │  │  │  ├── __init__.py
   │  │  │  ├── plate01.py
   │  │  │  └── runall.py
   │  │  ├── runall.py
   │  │  ├── solids
   │  │  │  ├── __init__.py
   │  │  │  └── runall.py
   │  │  └── trusses
   │  │     ├── __init__.py
   │  │     ├── runall.py
   │  │     ├── truss01.py
   │  │     ├── truss02.py
   │  │     ├── truss03.py
   │  │     └── truss04.py
   │  ├── materials
   │  │  ├── __init__.py
   │  │  ├── ElasticSection.py
   │  │  ├── FiberMaterial.py
   │  │  ├── Material.py
   │  │  ├── PlaneStrain.py
   │  │  ├── PlaneStress.py
   │  │  ├── SectionMaterial.py
   │  │  └── VonMises.py
   │  ├── plotter
   │  │  ├── __init__.py
   │  │  ├── AbstractPlotter.py
   │  │  ├── ElementPlotter.py
   │  │  ├── ElementPlotter3D.py
   │  │  ├── Plotter.py
   │  │  └── Plotter3D.py
   │  └── solver
   │     ├── __init__.py
   │     ├── LinearSolver.py
   │     ├── NewtonRaphsonSolver.py
   │     └── Solver.py
   ├── main.py
   └── WIP
      └── PlaneStress.py
