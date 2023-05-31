"""
==========================================================
A first 1D spring system
==========================================================

.. code::

                o---/\/\/\---o---/\/\/\---o
   x---/\/\/\---o                         o--->
                o-------/\/\/\/\/\--------o

   x ... fixed node
   o ... free node (moving in x-direction)


.. note::

    **This example uses**

    * :ref:`System class`
    * :ref:`Spring class`
    * :ref:`Node class`

"""

# %%
# Initialization
# ----------------
#

from femedu.domain import *
from femedu.elements.linear import Spring

# %%
# Building the model
# -------------------
# 1. Initializing a model

model = System()

# %%
# 2. Defining nodes
nd1 = Node(0.0, 0.0)
nd2 = Node(2.0, 0.0)
nd3 = Node(4.0, 0.0)
nd4 = Node(6.0, 0.0)

model.addNode(nd1, nd2, nd3, nd4)

# %%
# 3. Creating the springs
springA = Spring(nd1, nd2, 15)
springB = Spring(nd2, nd3, 10)
springC = Spring(nd3, nd4, 10)
springD = Spring(nd2, nd4, 10)

model.addElement(springA,springB,springC,springD)

# %%
# 4. Applying the essential boundary conditions
nd1.fixDOF('ux')

# %%
# 5. Applying loads
nd4.setLoad([2.0],['ux'])

# %%
# You may check your model any time by executing
model.report()

# %%
# Performing the analysis
# -----------------------
# 6. Assembly and solve
model.solve()

# %%
# 7. Check out displacements and internal forces
model.report()

# %%
# We can also create a force plot, though it doesn't look all that nice in 1D
model.beamValuePlot('f')

