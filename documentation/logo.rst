The Ocellaris logo
==================

The Ocellaris logo is based on a `picture of an Ocellaris clownfish 
<http://commons.wikimedia.org/wiki/File:Amphiprion_ocellaris_%28Clown_anemonefish%29_Nemo.jpg>`_
by Nick Hobgood from `Wikipedia 
<http://en.wikipedia.org/wiki/Ocellaris_clownfish>`_
licensed under `CC BY-SA
<http://creativecommons.org/licenses/by-sa/3.0/deed.en>`_.

The triangulation of the fish was made using an `online image triangulation tool
<http://snorpey.github.io/triangulation/>`_ by Georg Fischer. The image was first
prepared in the GIMP image editor; the background was replaced with white and the
fish was blurred and defeatured to look better when triangulated. Small features
made the resulting triangulation quite noisy so this step was necessary to get a
good result. The resulting image was then processed through the online tool and
exported from the web site to SVG format.

.. Triangulation tool settings 0 97 6 9

The SVG file from the web site was post-processed in Python with ElementTree to
remove all white triangles and then edited in Inkscape to adjust the color and
shape of a few triangles. The final logo image was then exported to PNG format
from Inkscape.

The computational mesh
----------------------

A computational mesh was made of the domain around the Ocellaris clownfish.
The PNG image file containing the logo was opened in the GIMP image editor
and the outline of the Ocellaris clownfish was selected and converted to a path.
This path was saved to SVG format in order to get the coordinates of the outline.
These coordinates were fed into Gmsh by generating a GEO file in Python and a 
triangulation of the outside of the fish was created in Gmsh MSH format.

The image below has been created in Python. The triangulation around the fish
was read from the Gmsh MSH format output file and converted to a Matplotlib 
triangulation. The logo image was read in by Matplotlib and superinposed on
a plot of the triangulation. This was the only hard step: to get the PNG and 
the plot of the triangulated mesh to align. If you look closely you will see
that the mesh does not follow the PNG image perfectly. I think it is good enough
for an evening hack.

.. image:: https://trlandet.bitbucket.io/figures/ocellaris_mesh_521.png
    :alt: Mesh with Ocellaris fish logo image superimposed 

Why?
----

Just for fun.

Both the triangulated look of the logo and the creation of the mesh for calculating
flow around the Ocellaris clownfish are references to the FEniCS logo and the flow
around a dolphin FEniCS demo that is included in the Dolfin FEniCS package.

Why this page?
..............

This page is mainly here to satisfy the CC-BY-SA license of the original image:

    **Attribution** — You must give appropriate credit, provide a link to the license,
    and indicate if changes were made. You may do so in any reasonable manner, but not
    in any way that suggests the licensor endorses you or your use.

    **ShareAlike** — If you remix, transform, or build upon the material, you must
    distribute your contributions under the same license as the original.

    No additional restrictions — You may not apply legal terms or technological
    measures that legally restrict others from doing anything the license permits.
    
The logo image is hence under the CC-BY-SA license.
