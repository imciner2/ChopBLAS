==============
:func:`chasum`
==============
For the vector :math:`x`, compute the quantity :math:`x_{out} = \sum_{i} \lvert x_{i} \rvert` with operation-level rounding.

-----------
Description
-----------



-----
Usage
-----

.. mat:function:: xout = chasum(x, ...)

   Compute :math:`x_{out} = \sum_{i} \lvert x_{i} \rvert` using the :func:`chop` rounding function and the global rounding options.

   :param x: The vector to compute over.
   :return: Computed value

.. mat:function:: xout = chasum(x, opts, ...)

   Compute :math:`x_{out} = \sum_{i} \lvert x_{i} \rvert` using the :func:`chop` rounding function and the rounding options given in `opts`.

   :param x: The vector to compute over.
   :param opts: Options to pass to :func:`chop`.
   :return: Computed value


Optional parameters
+++++++++++++++++++

--------
Examples
--------
